#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class SamInput : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            db.push_path( "sample_files", child.second );
        }
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );
        std::vector< std::string > genome_fastas( db.require_genome( db ));
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > sam_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Sam_ihandler_impl< IoHandlerIfstream > sam_reader( sam_paths );

        monitor.set_monitor( "Component SamInput", sam_paths.size() +2 );
        monitor.log( "Component SamInput", "Start" );

        std::mutex sam_mutex;
        ParaThreadPool sam_parallel_pool( sam_paths.size() );

        for( auto& sam_path : sam_paths )
        {
            std::string sample_name( get_sample_name( sam_path ));
            db.sam_samples.emplace_back( sample_name, std::vector< Sam<> >{} );
        }

        for( size_t smp = 0; smp < sam_paths.size(); ++smp )
        {
            sam_parallel_pool.job_post( [ smp, &db, &sam_paths, &sam_reader, &sam_mutex, &monitor, this ] ()
            {
                std::vector< Sam<> > sams;
                Sam<> sam;

			    while( true )
			    {
				    sam = sam_reader.get_next_entry( smp );

                    if( sam.eof_flag )
                    {
                        break;
                    }

                    sams.push_back( sam );
                }

                {
                    std::lock_guard< std::mutex > sam_lock( sam_mutex );
                    db.sam_samples[ smp ].second = std::move( sams );
                    monitor.log( "Component SamInput", ( db.sam_samples[ smp ].first ).c_str() );
                }
            });

        }

        sam_parallel_pool.flush_pool();
        monitor.log( "Component SamInput", "Complete" );
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
};

} // end of namespace component
} // end of namespace ago
