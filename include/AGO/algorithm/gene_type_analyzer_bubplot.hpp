#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBubplot
{
  public:

    GeneTypeAnalyzerBubplot()
    {}

    static void output_bubplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples
            )
    {
        for( auto& smp : bed_samples )
        {
            boost::filesystem::create_symlink(( "../LenDist/" + smp.first + ".tsv" ).c_str(), ( output_name + smp.first + ".tsv" ).c_str() );
        }
    }

    static void output_anno_sequence()
    {
    }

};

} // end of namespace algorithm
} // end of namespace ago
