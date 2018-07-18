#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <AGO/format/md_sam.hpp>
#include <mutex>

namespace ago {
namespace component {

class SamInput : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    bool allow_t2c;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );
        allow_t2c = p.get_optional< bool >( "allow_t2c" ).value_or( false );

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
        db.require_genome( db );
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > sam_paths( get_path_list_string( db.get_path_list( "sample_files" )));

        monitor.set_monitor( "Component SamInput", sam_paths.size() +2 );
        monitor.log( "Component SamInput", "Start" );

        std::mutex sam_mutex;
        ParaThreadPool sam_parallel_pool( sam_paths.size() );

        for( auto& sam_path : sam_paths )
        {
            std::string sample_name( get_sample_name( sam_path ));
            db.sam_samples.emplace_back( sample_name, std::vector< ago::format::MDSam<> >{} );
        }

        for( size_t smp = 0; smp < sam_paths.size(); ++smp )
        {
            sam_parallel_pool.job_post( [ smp, &db, &sam_paths, &sam_mutex, &monitor, this ] ()
            {
                std::ifstream input( sam_paths[ smp ] );
                std::vector< ago::format::MDSam<> > sams;
                std::string samline;

                while( std::getline( input, samline ))
			    {
                    if( samline.at(0) == '@' ) continue;
                    sams.push_back( ago::format::MDSam<>( samline, allow_t2c ));
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
