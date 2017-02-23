#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>
#include <boost/archive/text_oarchive.hpp>

namespace ago {
namespace component {

class Annotator : public engine::NamedComponent
{
    struct HitHandler
    {
        template< class DB_BED, class ANN_BED >
        static void run( DB_BED&& db_bed, ANN_BED&& ann_bed, int& db_idx )
        {
	        tuple2vector<
                std::vector<std::string>
                , 4
                , std::tuple_size< decltype(db_bed.data) >::value
            > ( db_bed.data , ann_bed.annotation_info_[db_idx] );

            if( std::get<4>( db_bed.data ) == "miRNA" )
            {
                for( auto& info : ann_bed.annotation_info_ )
                {
                    for( int i = 0; i < info.size(); i+=2 )
                    {
                        if( info[i] == "miRNA" && info[i+1].at( info[i+1].size()-1 ) != 'p' )
                        {
                            size_t db_mid  = std::get<1>( db_bed.data ) + (( std::get<2>( db_bed.data ) - std::get<1>( db_bed.data )) /2 );
                            size_t ann_mid = ann_bed.start_ + (( ann_bed.end_ - ann_bed.start_ ) /2 );

                            if( ann_mid > db_mid )
                            {
                                switch( ann_bed.strand_ )
                                {
                                    case '+' : info[i+1] = info[i+1] + "-3p"; break;
                                    case '-' : info[i+1] = info[i+1] + "-5p"; break;
                                }
                            }
                            else
                            {
                                switch( ann_bed.strand_ )
                                {
                                    case '+' : info[i+1] = info[i+1] + "-5p"; break;
                                    case '-' : info[i+1] = info[i+1] + "-3p"; break;
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    using Base = engine::NamedComponent;

    using BedFileReaderImpl =  FileReader_impl<
          Bed
        , std::tuple< std::string, uint32_t, uint32_t, char, std::string, std::string >
        , SOURCE_TYPE::IFSTREAM_TYPE
    >;

    using AnnotationTrait = Annotation<
          BedFileReaderImpl
        , AnnoIgnoreStrand::NO_IGNORE
        , AnnoType::INTERSET
        , HitHandler
    >;

    using Annotations = AnnotationSet <
          AnnotationRawBed<>
        , AnnotationTrait
    >;

    std::vector< std::string > annotation_files_;

    bool output_archive_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "annotation_files" ))
        {
            db.push_path( "annotation_files", child.second );
        }

        output_archive_ = p.get_optional< bool >( "output_archive" ).value_or( false );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );

        if( db.exist_path_tag( "annotation_files" ))
        {
            for( auto& bed : db.get_path_list( "annotation_files" ))
            {
                annotation_files_.emplace_back( bed.string() );
            }
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component Annotator", db.bed_samples.size() +3 );
        monitor.log( "Component Annotator", "Loading Annotation" );

        monitor.set_monitor( "	Loading Annotation", 2 );
        monitor.log( "	Loading Annotation", " ... " );

        Annotations annotator( annotation_files_ );
        monitor.log( "	Loading Annotation", " ... Done" );

        monitor.log( "Component Annotator", "Annotating Bed" );

        std::vector< std::ofstream > archive_outputs;

        if( output_archive_ )
        {
            for( size_t smp = 0; smp < db.bed_samples.size(); ++smp )
            {
                archive_outputs.push_back( std::move( std::ofstream(
                    db.output_dir().string() + db.bed_samples[ smp ].first + ".arc"
                )));
            }
        }

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.bed_samples.size() );

        std::pair< size_t, std::vector< AnnotationRawBed<> >> sample_bed_pair;

        for( size_t smp = 0; smp < db.bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [
                sample_bed_pair{ std::make_pair( smp, std::move( db.bed_samples[ smp ].second ))},
                &db, &archive_outputs, &annotator, &monitor, &smp_mutex, this ] () mutable
            {
                for( auto& anno_rawbed : sample_bed_pair.second )
                {
                    annotator.AnnotateAll( anno_rawbed );
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );
                    monitor.log( "Component Annotator", ( db.bed_samples[ sample_bed_pair.first ].first ).c_str() );

                    if( output_archive_ )
                    {
                        boost::archive::binary_oarchive archive_out( archive_outputs[ sample_bed_pair.first ] );
                        archive_out & sample_bed_pair.second;
                        archive_outputs[ sample_bed_pair.first ].close();
                    }

                    db.bed_samples[ sample_bed_pair.first ].second = std::move( sample_bed_pair.second );
                }
            });
        }

        smp_parallel_pool.flush_pool();
        Annotations::clear_database();

        monitor.log( "Component Annotator", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago