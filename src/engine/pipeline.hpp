#pragma once
#include <vector>
#include <AGO/engine/component.hpp>

namespace ago {
namespace engine {

template<class LIST>
class PipelineImpl : public LIST
{
protected:
    using Base = LIST;

public :
    PipelineImpl ( LIST&& component_list )
    : Base ( std::move( component_list ) )
    {}

    void operator() ()
    {
        for( auto& component_ptr : (*this) )
        {
            component_ptr->initialize();
            (*component_ptr)();
        }
    }
};

template<class LIST>
using Pipeline = PipelineImpl<LIST>;

template<class LIST>
auto make_pipeline( LIST&& component_list )
{
    static_assert ( std::is_rvalue_reference<LIST&&>::value, "");
    return Pipeline<std::decay_t<LIST>> ( 
        std::move( component_list ) 
    );
}

} // engine
} // ago
