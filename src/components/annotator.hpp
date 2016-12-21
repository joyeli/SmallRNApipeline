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

        monitor.set_monitor( "Component Annotator", 4 );
        monitor.log( "Component Annotator", "Start" );

        std::string output_dir = db.output_dir().string();

        monitor.set_monitor( "	Loading Annotation", 2 );
        monitor.log( "	Loading Annotation", " ... " );
        monitor.log( "Component Annotator", "Loading Annotation" );

        Annotations annotator( annotation_files_ );

        monitor.log( "	Loading Annotation", " ... Done" );
        monitor.set_monitor( "	Annotating Bed", db.bed_samples.size() +2 );
        monitor.log( "	Annotating Bed", "Start" );
        monitor.log( "Component Annotator", "Annotating Bed" );

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.bed_samples.size() );

        std::map< std::string, std::vector< AnnotationRawBed<> >> bed_samples_map;

        for( auto& sample : db.bed_samples )
        {
            smp_parallel_pool.job_post( [ sample, output_dir, &annotator, &monitor, &smp_mutex, &bed_samples_map, this ] () mutable
            {
                for( auto& anno_rawbed : sample.second )
                {
                    annotator.AnnotateAll( anno_rawbed );
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );

                    bed_samples_map.emplace( sample ); 
                    monitor.log( "	Annotating Bed", ( sample.first ).c_str() );

                    if( output_archive_ )
                    {
                        std::ofstream archive_output( output_dir + "/" + sample.first + ".arc" );
                        boost::archive::binary_oarchive archive_out( archive_output );

                        archive_out & sample.second;
                        archive_output.close();
                    }
                }
            });
        }

        smp_parallel_pool.flush_pool();
        Annotations::clear_database();

        for( auto& sample : db.bed_samples )
        {
            auto it = bed_samples_map.find( sample.first );
            sample.second.swap( it->second );
        }

        monitor.log( "	Annotating Bed", " ... Done" );
        monitor.log( "Component Annotator", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
