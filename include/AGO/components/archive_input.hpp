#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class ArchiveInput : public engine::NamedComponent
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

        std::vector< std::string > archive_paths( get_path_list_string( db.get_path_list( "sample_files" )));

        monitor.set_monitor( "Component ArchiveInput", archive_paths.size() +2 );
        monitor.log( "Component ArchiveInput", "Start" );

        for( auto& archive_path : archive_paths )
        {
            std::string sample_name( get_sample_name( archive_path ));
            db.bed_samples.emplace_back( sample_name, std::vector< AnnotationRawBed<> >{} );
        }

        for( size_t smp = 0; smp < archive_paths.size(); ++smp )
        {
            std::vector< AnnotationRawBed<> > annotation_rawbeds;

            std::ifstream archive( archive_paths[ smp ] );
            boost::archive::binary_iarchive archive_in( archive );

            archive_in & annotation_rawbeds;
            archive.close();

            db.bed_samples[ smp ].second = std::move( annotation_rawbeds );
            monitor.log( "Component ArchiveInput", ( db.bed_samples[ smp ].first ).c_str() );
        }

        monitor.log( "Component ArchiveInput", "Complete" );
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
