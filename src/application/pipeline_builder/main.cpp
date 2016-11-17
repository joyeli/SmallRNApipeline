#include <AGO/application/pipeline_builder/main.hpp>
#include <CPT/exception.hpp>
#include <CPT/exception/std_exception_json_dump.hpp>

int main ( int argc, const char* argv[] )
{
    ago::application::pipeline_builder::OptionParser option_parser( argc, argv );

    auto&& pipeline_builder( ago::application::pipeline_builder::make( option_parser ) );

    pipeline_builder();

    return 0;
}
