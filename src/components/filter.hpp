#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/annotator/filter.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>
#include <boost/archive/text_iarchive.hpp>

namespace ago {
namespace component {

class Filter : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    using FilterTypeList = boost::mpl::vector<
          boost::mpl::vector< boost::mpl::string< 'prot', 'ein_' >, boost::mpl::int_< 0 >, boost::mpl::char_< '^' >>
        , boost::mpl::vector< boost::mpl::string< 'misc', '_RNA' >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'nc'  , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 's'   , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'sn'  , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'sno' , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'sca' , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'ribo', 'zyme' >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'TR'  , '_'    >, boost::mpl::int_< 1 >, boost::mpl::char_< '^' >>
        , boost::mpl::vector< boost::mpl::string< 'IG'  , '_'    >, boost::mpl::int_< 1 >, boost::mpl::char_< '^' >>
        , boost::mpl::vector< boost::mpl::string< 'TE'  , 'C'    >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'r'   , 'msk'  >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>
        , boost::mpl::vector< boost::mpl::string< 'nc'  , 'rna'  >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 'sen' , 'se'   >, boost::mpl::int_< 1 >, boost::mpl::char_< '^' >>
        , boost::mpl::vector< boost::mpl::string< 'sen' , 'se'   >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 'eudo', 'gene' >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 'ansc', 'ript' >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 'r'   , 'RNA'  >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 't'   , 'RNA'  >, boost::mpl::int_< 1 >, boost::mpl::char_< '$' >>
        , boost::mpl::vector< boost::mpl::string< 'mi'  , 'RNA'  >, boost::mpl::int_< 0 >, boost::mpl::char_< '=' >>
    >;

    using Filters = FilterWorker< AnnotationRawBed<>, FilterTypeList >;

    bool archive_input_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

         archive_input_ = p.get_optional< bool >( "archive_input" ).value_or( false );

        if( archive_input_ )
        {
            for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
            {
                db.push_path( "sample_files", child.second );
            }
        }
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        if( archive_input_ )
        {
            auto& db( this->mut_data_pool() );

            std::vector< std::string > genome_fastas( db.require_genome( db ));
            std::vector< std::string > archive_paths( get_path_list_string( db.get_path_list( "sample_files" )));

            for( auto& archive_path : archive_paths )
            {
                std::string sample_name( get_sample_name( archive_path ));
                std::vector< AnnotationRawBed<> > annotation_rawbeds;

                std::ifstream archive( archive_path );
                boost::archive::binary_iarchive archive_in( archive );

                archive_in & annotation_rawbeds;
                archive.close();

                db.rawbed_samples.emplace_back( sample_name, annotation_rawbeds );
            }
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component Filter", 1 );

		Filters run_filter;

        monitor.set_monitor( "Filtering", db.rawbed_samples.size()+1 );

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.rawbed_samples.size() );
        std::map< std::string, std::vector< AnnotationRawBed<> >> rawbed_samples_map;

        for( auto& sample : db.rawbed_samples )
        {
            smp_parallel_pool.job_post( [ sample, &run_filter, &monitor, &smp_mutex, &rawbed_samples_map ] () mutable
            {
                for( auto& anno_rawbed : sample.second )
                {
                    anno_rawbed = run_filter.Filter( anno_rawbed );
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );
                    rawbed_samples_map.emplace( sample ); 

                    monitor.log( "Filtering", " ... " + sample.first );
                }
            });
        }

        smp_parallel_pool.flush_pool();

        for( auto& sample : db.rawbed_samples )
        {
            sample.second = rawbed_samples_map[ sample.first ];
        }

        monitor.log( "Filtering", " ... Complete" );
        monitor.log( "Component Filter", "Complete!!" );
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
