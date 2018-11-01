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

    void annobed_outputing( auto& output, auto& genome, auto& sample_beds, const auto& filter_ppm )
    {
        output << "Chr\tStart\tEnd\tStrand\tAlignCounts\tRawCounts\tReadCounts\tPPM\tRMSK\tLength\tTailLen\tSeq\tTail\tMM\tT2C\tType\tAnnoSeedMD\n";
    
        for( auto& anno : sample_beds )
        {
            if( anno.ppm_ < filter_ppm ) continue;
            // for( std::size_t i = 0; i < anno.annotation_info_.size(); ++i )
            {
                std::size_t i = 0; // do first priority
                if( i < anno.annotation_info_.size() )
                {
                    for( int j = 0; j < anno.annotation_info_[i].size(); j+=2 )
                    {
                        output
                            << anno.chromosome_ << "\t"
                            << anno.start_ << "\t"
                            << anno.end_ << "\t"
                            << anno.strand_ << "\t"
                            << anno.multiple_alignment_site_count_ << "\t"
                            << anno.reads_count_ << "\t"
                            << ((double)(anno.reads_count_) / (double)(anno.multiple_alignment_site_count_)) << "\t"
                            << anno.ppm_ << "\t"
                            << ( anno.is_filtered_ == 0 ? "N\t" : "Y\t" )
                            << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                            << (int)anno.tail_length_ << "\t"
                            << GeneTypeAnalyzerCounting::seqT2U( anno.getReadSeq( genome )) << "\t"
                            << ( GeneTypeAnalyzerCounting::seqT2U( anno.getTail() ) != "" ? GeneTypeAnalyzerCounting::seqT2U( anno.getTail() ) : "." ) << "\t"
                            << ( anno.md_map.size() != 0 ? GeneTypeAnalyzerCounting::seqT2U( anno.getMD() ): "." ) << "\t"
                            << ( anno.tc_set.size() != 0 ? GeneTypeAnalyzerCounting::seqT2U( anno.getTC() ): "." ) << "\t"
                            << anno.annotation_info_[i][j] << "\t"
                            << anno.annotation_info_[i][ j+1 ] << "_"
                            << GeneTypeAnalyzerCounting::seqT2U( anno.getReadSeq( genome ).substr(1,7) )
                            << ( GeneTypeAnalyzerCounting::seqT2U( anno.seed_md_tag ) != "" ? ( "|" + GeneTypeAnalyzerCounting::seqT2U( anno.seed_md_tag )) : "" )
                            << "\n"
                            ;
                    }
                }
            }
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
