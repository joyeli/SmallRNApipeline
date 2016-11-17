#pragma once
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <AGO/engine/components/named_component.hpp>
#include <CPT/logger.hpp>

namespace ago {
namespace component {

class FastqIn : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );

        db.fastq_min_length = p.get_optional< size_t >( "fastq_min_length" ).value_or( 15 );
        db.fastq_max_length = p.get_optional< size_t >( "fastq_max_length" ).value_or( 30 );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );
        db.start_shared_object_management();
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );

        if( db.fastq_samples.empty() )
        {
            for( auto& path : db.get_path_list( "sample_files" ))
            {
                db.require_fastq( path.string(), db );
            }
        }
    }
};

} // end of namespace component
} // end of namespace ago
