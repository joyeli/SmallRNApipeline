#pragma once
#include <AGO/engine/components/named_component.hpp>
// #include <pokemon/format/annotation_raw_bed.hpp>
#include <AGO/format/md_rawbed.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>
#include <boost/archive/text_oarchive.hpp>

namespace ago {
namespace component {

void output_annobed( auto& output, auto& genome, auto& sample_beds )
{
    output << "Chr\tStart\tEnd\tStrand\tAlignCounts\tRawCounts\tReadCounts\tPPM\tLength\tTailLen\tSeq\tTail\tMM\tT2C\tType\tAnnoSeedMD\n";

    for( auto& anno : sample_beds )
    {

        for( auto& info : anno.annotation_info_ )
        {
            if( info.size() == 0 )
            {
                output
                    << anno.chromosome_ << "\t"
                    << anno.start_ << "\t"
                    << anno.end_ << "\t"
                    << anno.strand_ << "\t"
                    << anno.multiple_alignment_site_count_ << "\t"
                    << anno.reads_count_ << "\t"
                    << ((double)(anno.reads_count_) / (double)(anno.multiple_alignment_site_count_)) << "\t"
                    << anno.ppm_ << "\t"
                    << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                    << (int)anno.tail_length_ << "\t"
                    << anno.getReadSeq( genome ) << "\t"
                    << ( anno.getTail() != "" ? anno.getTail() : "." ) << "\t"
                    << ( anno.md_map.size() != 0 ? anno.getMD() : "." ) << "\t"
                    << ( anno.tc_set.size() != 0 ? anno.getTC() : "." ) << "\t"
                    << ".\t"
                    << ".\n"
                    ;
            }
            else
            {
                for( int i = 0; i < info.size(); i+=2 )
                {
                    output
                        << anno.chromosome_ << "\t"
                        << anno.start_ << "\t"
                        << anno.end_ << "\t"
                        << anno.strand_ << "\t"
                        << anno.multiple_alignment_site_count_ << "\t"
                        << anno.reads_count_ << "\t"
                        << ((double)(anno.reads_count_) / (double)(anno.multiple_alignment_site_count_)) << "\t"
                        << anno.ppm_ << "\t"
                        << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                        << (int)anno.tail_length_ << "\t"
                        << anno.getReadSeq( genome ) << "\t"
                        << ( anno.getTail() != "" ? anno.getTail() : "." ) << "\t"
                        << ( anno.md_map.size() != 0 ? anno.getMD() : "." ) << "\t"
                        << ( anno.tc_set.size() != 0 ? anno.getTC() : "." ) << "\t"
                        << info[i] << "\t"
                        << info[ i+1 ] << "_"
                        << anno.getReadSeq( genome ).substr(1,7)
                        << ( anno.seed_md_tag != "" ? ( "|" + anno.seed_md_tag ) : "" )
                        << "\n"
                        ;
                }
            }
        }
    }
}

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

            if( std::get<4>( db_bed.data ).substr( 0, 5 ) == "miRNA" || std::get<4>( db_bed.data ) == "mirtron" )
            {
                for( auto& info : ann_bed.annotation_info_ )
                {
                    for( int i = 0; i < info.size(); i+=2 )
                    {
                        if(( info[i].substr( 0, 5 ) == "miRNA" || info[i] == "mirtron" )&& info[i+1].at( info[i+1].size()-1 ) != 'p' )
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
          // AnnotationRawBed<>
          ago::format::MDRawBed
        , AnnotationTrait
    >;

    std::vector< std::string > annotation_files_;

    bool output_archive_;
    bool output_annobed_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "annotation_files" ))
        {
            db.push_path( "annotation_files", child.second );
        }

        output_archive_ = p.get_optional< bool >( "output_archive" ).value_or( true );
        output_annobed_ = p.get_optional< bool >( "output_annobed" ).value_or( true );
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
        monitor.log( "	Loading Annotation", "..." );

        Annotations annotator( annotation_files_ );
        monitor.log( "	Loading Annotation", "Done" );

        monitor.log( "Component Annotator", "Annotating Bed" );

        std::vector< std::ofstream > archive_outputs;
        std::vector< std::ofstream > annobed_outputs;

        if( output_archive_ )
        {
            for( size_t smp = 0; smp < db.bed_samples.size(); ++smp )
            {
                archive_outputs.push_back( std::move( std::ofstream(
                    db.output_dir().string() + db.bed_samples[ smp ].first + ".arc"
                )));
            }
        }

        if( output_annobed_ )
        {
            for( size_t smp = 0; smp < db.bed_samples.size(); ++smp )
            {
                annobed_outputs.push_back( std::move( std::ofstream(
                    db.output_dir().string() + db.bed_samples[ smp ].first + "_annobed.text"
                )));
            }
        }

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.bed_samples.size() );

        // std::pair< size_t, std::vector< AnnotationRawBed<> >> sample_bed_pair;
        std::pair< size_t, std::vector< ago::format::MDRawBed >> sample_bed_pair;

        for( size_t smp = 0; smp < db.bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [
                sample_bed_pair{ std::make_pair( smp, std::move( db.bed_samples[ smp ].second ))},
                &db, &archive_outputs, &annobed_outputs, &annotator, &monitor, &smp_mutex, this ] () mutable
            {
                double sum = 0.0;

                for( auto& anno_rawbed : sample_bed_pair.second )
                {
                    annotator.AnnotateAll( anno_rawbed );
                    anno_rawbed.ppm_ = (double)(anno_rawbed.reads_count_) / (double)(anno_rawbed.multiple_alignment_site_count_);
                    sum += anno_rawbed.ppm_;
                }

                double ppm_rate = 1000000.0 / sum;

                for( auto& anno_rawbed : sample_bed_pair.second )
                {
                    anno_rawbed.ppm_ = anno_rawbed.ppm_ * ppm_rate;
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

                    if( output_annobed_ )
                    {
                        output_annobed( annobed_outputs[ sample_bed_pair.first ], db.genome_table, sample_bed_pair.second );
                        annobed_outputs[ sample_bed_pair.first ].close();
                    }

                    db.bed_samples[ sample_bed_pair.first ].second = std::move( sample_bed_pair.second );
                }
            });
        }


        smp_parallel_pool.flush_pool();
        Annotations::clear_database();

        output_statistic( db.output_dir().string(), db.bed_samples );
        monitor.log( "Component Annotator", "Complete" );
    }

    void output_statistic( const auto& output_dir, auto& bed_samples )
    {
        std::size_t anno_count = 0;
        std::set< std::string > bioindex;
        std::map< std::string, double > temp_map;

        std::vector< std::map< std::string, double >>
            statistic_samples( bed_samples.size(), std::map< std::string, double >() );

        std::vector< double > statistic_totals( bed_samples.size(), 0.0 );

        for( size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& anno : bed_samples[ smp ].second )
            {
                for( auto& info : anno.annotation_info_ )
                {
                    for( int i = 0; i < info.size(); i+=2 )
                    {
                        if( temp_map.find( info[i] ) == temp_map.end() )
                        {
                            bioindex.emplace( info[i] );
                            temp_map[ info[i] ] = 0.0;
                        }

                        temp_map[ info[i] ] += (double)(anno.reads_count_) / (double)(anno.multiple_alignment_site_count_);
                        anno_count += 1;
                    }
                }

                for( auto& biotype : temp_map )
                {
                    statistic_samples[ smp ][ biotype.first ] += biotype.second / (double)(anno_count);
                    statistic_totals[ smp ] += biotype.second / (double)(anno_count);
                }

                temp_map.clear();
                anno_count = 0;
            }
        }

        std::ofstream output( output_dir + "annotated.text" );

        output << "Sample";

        for( auto& smp : bed_samples )
        {
            output << "\t" << smp.first;
        }

        output << "\n";

        for( auto& biotype : bioindex )
        {
            output << biotype;

            for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                output << std::fixed << std::setprecision( 0 ) << "\t" << statistic_samples[ smp ][ biotype ];
            }
            
            output << "\n";
        }

        output << "Total";

        for( size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
                output << std::fixed << std::setprecision( 0 ) << "\t" << statistic_totals[ smp ];
        }

        output << "\n";
        output.close();
    }
};

} // end of namespace component
} // end of namespace ago
