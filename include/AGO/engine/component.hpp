#pragma once
#include <AGO/engine/data_pool.hpp>
#include <memory>
#include <boost/property_tree/ptree.hpp>

/**
 * the framework of component
 * do not add any logical code in this section
 */

namespace ago {
namespace engine {

class ComponentImpl
{
  protected: 

    const DataPool& data_pool_;
  public:

    ComponentImpl ( 
          const DataPool& data_pool
        , const bpt::ptree& schema_node
    )
    : data_pool_ ( data_pool )
    {
    }

    static  void build_rule () {}
    virtual void initialize() {}
    virtual void config ( const bpt::ptree& parameter ) {}

    virtual void start      () {}
    virtual void finish     () {}

    virtual void operator() ()
    {
        this->start();
        this->finish();
    }

    DataPool& mut_data_pool()
    {
        return const_cast<DataPool&>(data_pool_);
    }
    virtual ~ComponentImpl() {}
};

using Component = ComponentImpl;

} // engine
} // ago

#include <AGO/engine/components/named_component.hpp>

namespace ago {
namespace engine {

using ComponentPtr = std::unique_ptr<Component>;
using NamedComponent = components::NamedComponent;

} // engine
} // ago
