#pragma once
#include <AGO/engine/component.hpp>
#include <memory>
#include <string>
#include <AGO/engine/data_pool.hpp>
#include <AGO/components/test.hpp>
#include <AGO/components/fastq_in.hpp>
#include <AGO/components/aligner_tailor.hpp>

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
        ID_MAP_TYPE( "Test" , component::Test )
        ID_MAP_TYPE( "FastqIn" , component::FastqIn )
        ID_MAP_TYPE( "AlignerTailor" , component::AlignerTailor )
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
