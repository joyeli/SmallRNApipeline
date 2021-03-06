#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>
#include <AGO/format/md_rawbed.hpp>
#include <AGO/format/md_sam.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class Formats
{
    using VcfType = std::tuple< std::string, std::size_t, std::string, std::string, double, std::string >;
    //                              chr         pos         ref             alt     qul         id
  
  public:

    std::vector< std::pair<
        std::string,
        std::vector< Fastq<> >
    >> fastq_samples;

    std::vector< std::pair<
        std::string,
        std::vector< ago::format::MDSam<> >
    >> sam_samples;

    std::vector< std::pair<
        std::string,
        std::vector< ago::format::MDRawBed >
    >> bed_samples;

    template< class DB >
    Formats( DB& db )
    {
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
