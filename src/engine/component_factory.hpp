#pragma once
#include <AGO/engine/component.hpp>
#include <memory>
#include <string>
#include <AGO/engine/data_pool.hpp>
#include <AGO/components/tailor_fastq.hpp>
#include <AGO/components/annotator.hpp>
#include <AGO/components/filter.hpp>
#include <AGO/components/analyzer.hpp>
#include <AGO/components/meta_analyzer.hpp>
#include <AGO/components/visualization.hpp>

namespace ago {
namespace engine {

#define ID_MAP_TYPE( ID, TYPE ) \
    else if ( identifier == ID ) \
    { \
        auto tmp = ComponentPtr( \
            new TYPE(data_pool_, schema_node) \
        ); \
        tmp->config(schema_node); \
        return tmp; \
    }

class ComponentFactory
{
  public : 

    const DataPool& data_pool_;

    ComponentFactory( const DataPool& data_pool )
    : data_pool_( data_pool )
    {}

    ComponentPtr create_by_identifier ( 
          const std::string& identifier 
        , const bpt::ptree& schema_node
    )
    {

        if( identifier == "" ){ /* TODO a error handle */ }
        ID_MAP_TYPE( "TailorFastq" , component::TailorFastq )
        ID_MAP_TYPE( "Annotator" , component::Annotator )
        ID_MAP_TYPE( "Filter" , component::Filter )
        ID_MAP_TYPE( "Analyzer", component::Analyzer )
        ID_MAP_TYPE( "MetaAnalyzer", component::MetaAnalyzer )
        ID_MAP_TYPE( "Visualization", component::Visualization )
        /* TODO ADD COMPONENT HERE */
        else { /* TODO another error handle */ }

        return nullptr;
    }

    template<class T>
    auto operator () ( const T& component_schema )
    {
        auto name = component_schema
            .second
            .template get<std::string> ("name"); 
        return create_by_identifier ( 
              name
            , component_schema.second
        );
    }
};

#undef ID_MAP_TYPE

} // engine
} // ago
