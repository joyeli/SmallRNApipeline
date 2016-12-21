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

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component Filter", db.bed_samples.size() +2 );
        monitor.log( "Component Filter", "Start" );

		Filters run_filter;

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.bed_samples.size() );

        std::map< std::string, std::vector< AnnotationRawBed<> >> rawbed_samples_map;

        for( auto& sample : db.bed_samples )
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
                    monitor.log( "Component Filter", ( sample.first ).c_str() );
                }
            });
        }

        smp_parallel_pool.flush_pool();

        for( auto& sample : db.bed_samples )
        {
            sample.second = rawbed_samples_map[ sample.first ];
        }

        monitor.log( "Component Filter", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
