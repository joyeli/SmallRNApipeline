#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerAnnobed
{
  public:

    GeneTypeAnalyzerAnnobed()
    {}

    void annobed_outputing( auto& output, auto& genome, auto& sample_beds )
    {
        double ppm = GeneTypeAnalyzerCounting::get_ppm( sample_beds );
        output << "Chr\tStart\tEnd\tStrand\tAlignCounts\tRawCounts\tReadCounts\tPPM\tRMSK\tLength\tTailLen\tSeq\tTail\tType\tAnnoSeed\n";
    
        for( auto& anno : sample_beds )
        {
            for( auto& info : anno.annotation_info_ )
            {
                if( info.size() == 0 )
                {
                    output
                        << anno.chromosome_ << "\t"
                        << anno.start_ << "\t"
                        << anno.end_ << "\t"
                        << anno.strand_ << "\t"
                        << anno.multiple_alignment_site_count_ << "\t"
                        << anno.reads_count_ << "\t"
                        << anno.reads_count_       / anno.multiple_alignment_site_count_ << "\t"
                        << anno.reads_count_ * ppm / anno.multiple_alignment_site_count_ << "\t"
                        << ( anno.is_filtered_ == 0 ? "N\t" : "Y\t" )
                        << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                        << (int)anno.tail_length_ << "\t"
                        << anno.getReadSeq( genome ) << "\t"
                        << ( anno.getTail() != "" ? anno.getTail() : "." ) << "\t"
                        << ".\t"
                        << ".\n"
                        ;
                }
                else
                {
                    for( int i = 0; i < info.size(); i+=2 )
                    {
                        output
                            << anno.chromosome_ << "\t"
                            << anno.start_ << "\t"
                            << anno.end_ << "\t"
                            << anno.strand_ << "\t"
                            << anno.multiple_alignment_site_count_ << "\t"
                            << anno.reads_count_ << "\t"
                            << anno.reads_count_       / anno.multiple_alignment_site_count_ << "\t"
                            << anno.reads_count_ * ppm / anno.multiple_alignment_site_count_ << "\t"
                            << ( anno.is_filtered_ == 0 ? "N\t" : "Y\t" )
                            << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                            << (int)anno.tail_length_ << "\t"
                            << anno.getReadSeq( genome ) << "\t"
                            << ( anno.getTail() != "" ? anno.getTail() : "." ) << "\t"
                            << info[i] << "\t"
                            << info[ i+1 ] << "_"
                            << anno.getReadSeq( genome ).substr(1,7) << "\n"
                            ;
                    }
                }
            }
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
