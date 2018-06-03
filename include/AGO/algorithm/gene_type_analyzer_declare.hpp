#pragma once
#include <AGO/format/md_rawbed.hpp>

namespace ago {
namespace algorithm {

//                                   annotation                 length,  count
using CountingTableType = std::map< std::string, std::map< std::size_t, double >>;

//                                       annotation_index             length_index
using AnnoLengthIndexType = std::pair< std::set< std::string >, std::set< std::size_t >>;

// using BedSampleType = std::pair< std::string, std::vector< AnnotationRawBed<> >>;    
using BedSampleType = std::pair< std::string, std::vector< ago::format::MDRawBed >>; 

} // end of namespace algorithm
} // end of namespace ago
