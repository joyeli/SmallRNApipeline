#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <pokemon/converter/sam2rawbed.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class TailorFastqToBed : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t fastq_min_length_;
    size_t fastq_max_length_;

    size_t align_min_length_;
    size_t align_max_length_;

    size_t align_job_number_;
    size_t align_thread_num_;
    size_t align_limit_algn_;

    Aligner< Aligner_trait< Aligner_types::BWT_Aligner, ParallelTypes::M_T, void >> tailor_;

    bool output_sam_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            db.push_path( "sample_files", child.second );
        }

        std::string tailor_index = p.get_optional< std::string >( "tailor_index" ).value_or( "" );

        if( tailor_index != "" && !db.exist_path_tag( "tailor_index" ))
        {
            db.push_path( "tailor_index", tailor_index );
        }

        fastq_min_length_ = p.get_optional< size_t >( "fastq_min_length" ).value_or( 15 );
        fastq_max_length_ = p.get_optional< size_t >( "fastq_max_length" ).value_or( 30 );

        align_min_length_ = p.get_optional< size_t >( "align_min_length" ).value_or( 18 );
        align_max_length_ = p.get_optional< size_t >( "align_max_length" ).value_or( 30 );

        align_job_number_ = p.get_optional< size_t >( "align_job_number" ).value_or( 1000000 );
        align_thread_num_ = p.get_optional< size_t >( "align_thread_num" ).value_or( 20 );
        align_limit_algn_ = p.get_optional< size_t >( "align_limit_algn" ).value_or( 100 );

        output_sam_ = p.get_optional< bool >( "output_sam" ).value_or( false );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        if( !db.exist_path_tag( "tailor_index" ))
        {
            auto path = db.output_dir();
            path += "tailor";
            db.push_path( "tailor_index", path );
        }

        std::vector< std::string > genome_fastas( db.require_genome( db ));
        std::string tailor_index( db.require_tailor_index( genome_fastas, db ));

        if( db.is_tailor_index_build )
        {
            monitor.set_monitor( "Loading Index", 2 );
            monitor.log( "Loading Index", "Start" );

            tailor_.load_table( tailor_index );
            monitor.log( "Loading Index", "Complete" );
        }
        else
        {
            monitor.set_monitor( "Building Index", 3 );
            monitor.log( "Building Index", "Building ... " );

            tailor_.build( genome_fastas, tailor_index );
            monitor.log( "Building Index", "Loading ... " );

            tailor_.load_table( tailor_index );
            db.is_tailor_index_build = true;

            monitor.log( "Building Index", "Complete" );
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > fastq_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fastq_ihandler_impl< IoHandlerIfstream > fastq_reader( fastq_paths );

        monitor.set_monitor( "Component TailorFastqToBed", fastq_paths.size() +2 );
        monitor.log( "Component TailorFastqToBed", "Start" );

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( fastq_paths.size() );


        for( size_t id = 0; id < fastq_paths.size(); ++id )
        {
            smp_parallel_pool.job_post( [ id, &db, &fastq_paths, &fastq_reader, &smp_mutex, &monitor, this ] ()
            {
                std::string sample_name( get_sample_name( fastq_paths[ id ] ));

                std::mutex ali_mutex;
                ParaThreadPool ali_parallel_pool( align_thread_num_ );
                
                std::ofstream output_sam;

                if( output_sam_ )
                {
                    output_sam.open( db.output_dir().string() + "/" + sample_name + ".sam" );
                }

                std::vector< Sam<> > sams{};

                size_t thread_count = 0;
                bool break_flag = false;

                std::vector< Fastq<> > fastqs;

                while( true )
                {
                    for( size_t job = 0; job < align_job_number_; ++job )
                    {
                        Fastq<> fastq( fastq_reader.get_next_entry( id ));

                        if( fastq.eof_flag )
                        {
                            break_flag = true;
                            break;
                        }

                        if( fastq.getSeq().size() < fastq_min_length_ ||
                            fastq.getSeq().size() > fastq_max_length_ ||
                            n_check( fastq )
                          )
                            continue;

                        fastqs.push_back( fastq );
                    }

                    ali_parallel_pool.job_post( [ fastqs, &ali_mutex, &sams, &output_sam, this ] () 
                    {
                        std::vector< Sam<> > sams_para{};
                        std::map< int, std::vector< Fastq<> >> input;
                        std::vector< Sam<> >* sams_tmp;

                        for( auto& fastq: fastqs )
                        {
                            input.emplace( 0, std::vector< Fastq<> >{ fastq });
                            sams_tmp = tailor_.search( &input, 1, align_min_length_, align_limit_algn_, align_max_length_ );

                            if( !sams_tmp->empty() )
                            {
                                std::move( sams_tmp->begin(), sams_tmp->end(), std::back_inserter( sams_para ));
                            }

                            input.clear();
                            sams_tmp->clear();
                        }

                        {
                            std::lock_guard< std::mutex > ali_lock( ali_mutex );

                            if( output_sam_ )
                            {
                                for( auto& sam : sams_para )
                                {
                                    output_sam << sam;
                                }
                            }

                            std::move( sams_para.begin(), sams_para.end(), std::back_inserter( sams ));
                        }
                    });

                    thread_count++;
                    fastqs.clear();

                    if( thread_count >= align_thread_num_ )
                    {
                        thread_count = 0;
                        ali_parallel_pool.flush_pool();
                    }

                    if( break_flag )
                    {
                        ali_parallel_pool.flush_pool();

                        if( output_sam_ )
                        {
                            output_sam.close();
                        }

                        break;
                    }
                }

                for( auto sam = sams.begin(); sam != sams.end(); ++sam )
                {
                    std::get<3>( sam->data ) = std::get<3>( sam->data ) - 1;
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );

                    Sam2RawBed< std::vector< Sam<> >* > sam2bed;
                    auto raw_beds( sam2bed.Convert( &sams ));

                    std::vector< AnnotationRawBed<> > annotation_rawbeds;

                    for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
                    {
                        annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
                    }

                    sam2bed.rawbed_map_->clear();
                    sam2bed.rawbed_map2_->clear();

                    db.bed_samples.emplace_back( sample_name, annotation_rawbeds );
                    monitor.log( "Component TailorFastqToBed", ( sample_name ).c_str() );
                }
            });
        }

        smp_parallel_pool.flush_pool();
        monitor.log( "Component TailorFastqToBed", "Complete" );
    }

    std::vector< std::string > get_path_list_string( const std::vector< boost::filesystem::path >& paths )
    {
        std::vector< std::string > res;

        for( auto& path : paths )
        {
            res.emplace_back( path.string() );
        }

        return res;
    }

    std::string get_sample_name( const std::string& path )
    {
        std::vector< std::string > path_file;
        boost::iter_split( path_file, path, boost::algorithm::first_finder( "/" ));

        std::vector< std::string > sample;
        boost::iter_split( sample, path_file[ path_file.size()-1 ], boost::algorithm::first_finder( "." ));

        return sample[0];
    }

    bool n_check( const Fastq<>& fastq )
    {
        for( auto& base : fastq.getSeq() )
        {
            switch( base )
            {
                case 'A': break;
                case 'T': break;
                case 'C': break;
                case 'G': break;
                default : return true;
            }
        }

        return false;
    }
};

} // end of namespace component
} // end of namespace ago
