#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class Formats
{
  public:

    std::vector< std::pair<
        std::string,
        std::vector< Fastq<> >
    >> fastq_samples;

    std::vector< std::pair<
        std::string,
        std::vector< Sam<> >
    >> sam_samples;

    std::vector< std::pair<
        std::string,
        std::vector< AnnotationRawBed<> >
    >> bed_samples;

    std::vector< std::pair<
        std::string,
        std::vector< double >
    >> statistic_samples;

    template< class DB >
    Formats( DB& db )
    {
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
