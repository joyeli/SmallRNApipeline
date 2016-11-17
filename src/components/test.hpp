#pragma once
#include <AGO/engine/components/named_component.hpp>

namespace ago {
namespace component {

class Test : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
    }
};

} // end of namespace component
} // end of namespace ago
