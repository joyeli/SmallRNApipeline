#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <AGO/algorithm/quantile.hpp>
#include <CPT/logger.hpp>

namespace ago {
namespace engine {
namespace data_pool {


class MetaAnalyzer
{
  public:

    using anno_len_samples = std::map< std::string, std::map< std::string, std::vector< double >>>;
    //                                  SampleName              Anno                Length

    std::vector< std::pair< std::string, anno_len_samples >> quantile_result_samples;
    //                      AnalysisType

    template< class DB >
    MetaAnalyzer( DB& db )
    {
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
