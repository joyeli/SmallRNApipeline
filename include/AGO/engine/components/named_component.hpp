#pragma once 
#include <AGO/engine/component.hpp>
#include <CCD/utility/language.hpp>

namespace ago {
namespace engine {
namespace components {

class NamedComponentImpl : public ago::engine::Component
{
protected :
    virtual void config_parameters ( const bpt::ptree& parameters )
    {}
public :
    using Base = ago::engine::Component;
    NamedComponentImpl( 
          const DataPool& data_pool
        , const bpt::ptree& schema_node
    )
    :Base ( data_pool, schema_node )
    {
    }
    virtual void config ( const bpt::ptree& node_schema ) override
    {
        auto&& parameters = node_schema
            .get_child("parameter");
        this->config_parameters(parameters);
    }
};

CREATE_DERIVED_TYPE2( NamedComponent, NamedComponentImpl );

} // components
} // engine
} // ago
