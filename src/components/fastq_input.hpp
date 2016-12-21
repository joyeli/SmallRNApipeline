#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class FastqInput : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t fastq_min_length_;
    size_t fastq_max_length_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            db.push_path( "sample_files", child.second );
        }

        fastq_min_length_ = p.get_optional< size_t >( "fastq_min_length" ).value_or( 15 );
        fastq_max_length_ = p.get_optional< size_t >( "fastq_max_length" ).value_or( 30 );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > fastq_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fastq_ihandler_impl< IoHandlerIfstream > fastq_reader( fastq_paths );

        monitor.set_monitor( "Component FastqInput", fastq_paths.size()+2 );
        monitor.log( "Component FastqInput", "Start" );

        std::mutex fq_mutex;
        ParaThreadPool fq_parallel_pool( fastq_paths.size() );

        for( auto& fastq_path : fastq_paths )
        {
            std::string sample_name( get_sample_name( fastq_path ));
            db.fastq_samples.emplace_back( sample_name, std::vector< Fastq<> >{} );
        }

        for( size_t smp = 0; smp < fastq_paths.size(); ++smp )
        {
            fq_parallel_pool.job_post( [ smp, &db, &fastq_paths, &fastq_reader, &fq_mutex, &monitor, this ] ()
            {
                std::vector< Fastq<> > fastqs;
                Fastq<> fastq;

                while( true )
                {
                    fastq = fastq_reader.get_next_entry( smp );

                    if( fastq.eof_flag )
                    {
                        break;
                    }

                    if( fastq.getSeq().size() < fastq_min_length_ ||
                        fastq.getSeq().size() > fastq_max_length_ ||
                        n_check( fastq )
                    )
                    {
                        continue;
                    }

                    fastqs.push_back( fastq );
                }

                {
                    std::lock_guard< std::mutex > fq_lock( fq_mutex );
                    db.fastq_samples[ smp ].second = std::move( fastqs );
                    monitor.log( "Component FastqInput", ( db.fastq_samples[ smp ].first ).c_str() );
                }
            });
        }

        fq_parallel_pool.flush_pool();
        monitor.log( "Component FastqInput", "Complete" );
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
