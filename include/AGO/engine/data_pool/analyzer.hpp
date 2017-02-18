#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <CPT/logger.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class Analyzer
{
  public:

    //                              Each    LenMir_Map  Anno_Seed           All_Length      Length
    using AnalyzerResultType = std::vector< std::map< std::string, std::map< std::string, double >>>;

    AnalyzerResultType analyzer_result;
    std::vector< std::pair< std::string, AnalyzerResultType >> analyzer_result_samples;

    template< class DB >
    Analyzer( DB& db )
    {
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
