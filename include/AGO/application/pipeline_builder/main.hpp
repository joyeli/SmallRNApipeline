#pragma once
#include <CPT/config_reader/config_reader.hpp>
#include <CPT/option_parser.hpp>
#include <AGO/engine/data_pool.hpp>
#include <AGO/engine/component.hpp>
#include <AGO/engine/pipeline.hpp>
#include <AGO/engine/component_factory.hpp>
#include <CCD/utility/language.hpp>

namespace ago {
namespace application {
namespace pipeline_builder {

class OptionParser : public cpt::OptionParser
{
public:
    std::ifstream pipeline_schema_stream_;
    std::ofstream result_schema_stream_;

    OptionParser( int argc, char const * argv[] )
    {
        cpt::po::options_description desc( "Allowed options" );
        desc.add_options()
            ( "help,h", "show help message" )
            ( "input,i", cpt::po::value< std::string >(), "input json file" )
        ;

        cpt::po::store( cpt::po::parse_command_line( argc, argv, desc ), vm );
        cpt::po::notify( vm );

        if( vm.count( "help" ))
        {
            std::cout << desc << "\n";
            exit(1);
        }

        std::string input = "";
        get_parameter( "input", input );
        pipeline_schema_stream_.open( input == "" ? "/dev/stdin" : input );
    }
};

namespace bpt = boost::property_tree;
template< class OPTION_PARSER >
class Main
{
    OPTION_PARSER args_;
    ago::engine::DataPool data_pool_;

  public:

    Main( OPTION_PARSER&& args )
    : args_( std::forward<OPTION_PARSER>( args ))
    , data_pool_( args_ )
    {}

    DISABLE_COPY( Main );
    DEFAULT_MOVE( Main );

    void build_rule()
    {}

    auto build_pipeline()
    {
        std::vector< ago::engine::ComponentPtr > res;
        ago::engine::ComponentFactory component_factory( data_pool_ );

        auto&& pipeline_schema( data_pool_.pipeline_schema() ); 
        auto children( pipeline_schema.get_child( "pipeline" ));

        for( auto& child : children )
        {
            res.emplace_back( component_factory( child ));
        }

        return ago::engine::make_pipeline( std::move( res ));
    }

    void operator()()
    {
        auto pipeline( build_pipeline() );
        pipeline();
    }
};

template< class OPTION_PARSER >
Main< OPTION_PARSER > make( OPTION_PARSER&& option_parser )
{
    return Main< OPTION_PARSER >( std::forward< OPTION_PARSER >( option_parser ));
}

} // ago
} // application
} // pipeline_builder
