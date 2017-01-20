#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <mutex>

namespace ago {
namespace component {

class TailorAlign : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t align_min_length_;
    size_t align_max_length_;
    size_t align_limit_algn_;

    size_t thread_num_;

    Aligner< Aligner_trait< Aligner_types::BWT_Aligner, ParallelTypes::M_T, void >> tailor_;

    bool output_sam_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        std::string tailor_index = p.get_optional< std::string >( "tailor_index" ).value_or( "" );

        if( tailor_index != "" && !db.exist_path_tag( "tailor_index" ))
        {
            db.push_path( "tailor_index", tailor_index );
        }


        align_min_length_ = p.get_optional< size_t >( "align_min_length" ).value_or( 18 );
        align_max_length_ = p.get_optional< size_t >( "align_max_length" ).value_or( 30 );
        align_limit_algn_ = p.get_optional< size_t >( "align_limit_algn" ).value_or( 100 );

        thread_num_ = p.get_optional< size_t >( "thread_num" ).value_or( 32 );
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

        monitor.set_monitor( "Component TailorAlign", 4 );
        monitor.log( "Component TailorAlign", "Start" );

        std::vector< Fastq<> > fastqs;
        std::vector< std::pair< size_t, std::vector< Fastq<> >>> sample_fastqs_pool;

        monitor.log( "Component TailorAlign", "Loading Fastq" );
        monitor.set_monitor( "	Loading Fastq", db.fastq_samples.size() +2 );
        monitor.log( "	Loading Fastq", "Start" );

        size_t task_number = ( db.fastq_samples.size() / thread_num_ );

        for( size_t smp = 0; smp < db.fastq_samples.size(); ++smp )
        {
            for( auto& fastq : db.fastq_samples[ smp ].second )
            {
                fastqs.push_back( std::move( fastq ));

                if( fastqs.size() == task_number )
                {
                    sample_fastqs_pool.emplace_back( smp, std::move( fastqs ));
                    fastqs.clear();
                }
            }

            monitor.log( "	Loading Fastq", " ... " );
            sample_fastqs_pool.emplace_back( smp, std::move( fastqs ));
            fastqs.clear();
        }

        monitor.log( "	Loading Fastq", " ... Done" );
        monitor.log( "Component TailorAlign", "Aligning Fastq" );
        monitor.set_monitor( "	Aligning Fastq", sample_fastqs_pool.size() +2 );
        monitor.log( "	Aligning Fastq", "Start" );

        for( size_t smp = 0; smp < db.fastq_samples.size(); ++smp )
        {
            db.sam_samples.emplace_back( db.fastq_samples[ smp ].first, std::vector< Sam<> >{} );
        }

        std::vector< std::ofstream > sam_outputs;

        if( output_sam_ )
        {
            for( size_t smp = 0; smp < db.fastq_samples.size(); ++smp )
            {
                sam_outputs.push_back( std::move( std::ofstream(
                    db.output_dir().string() + "/" + db.fastq_samples[ smp ].first + ".sam"
                )));
            }
        }

        std::mutex ali_mutex;
        ParaThreadPool ali_parallel_pool( thread_num_ );
        size_t thread_count = 0;

        std::pair< size_t, std::vector< Fastq<> >> sample_fastqs_pair;

        for( size_t job = 0; job < sample_fastqs_pool.size(); ++job )
        {
            ali_parallel_pool.job_post( [
                sample_fastqs_pair{ std::move( sample_fastqs_pool[ job ] )},
                    &db, &sam_outputs, &ali_mutex, &monitor, this ] () mutable
            {
                std::vector< Sam<> > sams_para{};
                std::vector< Sam<> >* sams_tmp;

                std::map< int, std::vector< Fastq<> >> map_fastqs;
                map_fastqs.emplace( 0, std::vector< Fastq<> >{} );

                for( auto& fastq : sample_fastqs_pair.second )
                {
                    map_fastqs[0].push_back( std::move( fastq ));
                    sams_tmp = tailor_.search( &map_fastqs, 1, align_min_length_, align_limit_algn_, align_max_length_ );

                    if( !sams_tmp->empty() )
                    {
                        for( auto sam_it = sams_tmp->begin(); sam_it != sams_tmp->end(); ++sam_it )
                            std::get<3>( sam_it->data ) = std::get<3>( sam_it->data ) - 1;

                        std::move( sams_tmp->begin(), sams_tmp->end(), std::back_inserter( sams_para ));
                    }

                    map_fastqs[0].clear();
                    sams_tmp->clear();
                }

                {
                    std::lock_guard< std::mutex > ali_lock( ali_mutex );
                    monitor.log( "	Aligning Fastq", " ... " );

                    if( !sams_para.empty() )
                    {
                        if( output_sam_ )
                        {
                            for( auto& sam : sams_para )
                            {
                                sam_outputs[ sample_fastqs_pair.first ] << sam;
                            }
                        }

                        std::move(
                            sams_para.begin(),
                            sams_para.end(),
                            std::back_inserter( db.sam_samples[ sample_fastqs_pair.first ].second )
                        );
                    }
                }
            });

            thread_count++;

            if( thread_count >= thread_num_ )
            {
                thread_count = 0;
                ali_parallel_pool.flush_pool();
            }
        }

        ali_parallel_pool.flush_pool();

        monitor.log( "	Aligning Fastq", " ... Done" );

        if( output_sam_ )
        {
            for( auto& sam_output : sam_outputs )
            {
                sam_output.close();
            }
        }

        db.fastq_samples.clear();

        monitor.log( "Component TailorAlign", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
