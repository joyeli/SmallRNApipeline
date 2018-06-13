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
        output << "Chr\tStart\tEnd\tStrand\tAlignCounts\tRawCounts\tReadCounts\tPPM\tRMSK\tLength\tTailLen\tSeq\tTail\tMM\tT2C\tType\tAnnoSeedMD\n";
    
        for( auto& anno : sample_beds )
        {
            for( std::size_t i = 0; i < anno.annotation_info_.size(); ++i )
            {
                for( int j = 0; j < anno.annotation_info_[i].size(); j+=2 )
                {
                    if( anno.is_on_biotype_list_[i][j] )
                    {
                        output
                            << anno.chromosome_ << "\t"
                            << anno.start_ << "\t"
                            << anno.end_ << "\t"
                            << anno.strand_ << "\t"
                            << anno.multiple_alignment_site_count_ << "\t"
                            << anno.reads_count_ << "\t"
                            << (double)(anno.reads_count_) / (double)(anno.multiple_alignment_site_count_) << "\t"
                            << anno.reads_count_ * ppm / anno.multiple_alignment_site_count_ << "\t"
                            << ( anno.is_filtered_ == 0 ? "N\t" : "Y\t" )
                            << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                            << (int)anno.tail_length_ << "\t"
                            << anno.getReadSeq( genome ) << "\t"
                            << ( anno.getTail() != "" ? anno.getTail() : "." ) << "\t"
                            << ( anno.md_map.size() != 0 ? anno.getMD() : "." ) << "\t"
                            << ( anno.tc_set.size() != 0 ? anno.getTC() : "." ) << "\t"
                            << anno.annotation_info_[i][j] << "\t"
                            << anno.annotation_info_[i][ j+1 ] << "_"
                            << anno.getReadSeq( genome ).substr(1,7)
                            << ( anno.seed_md_tag != "" ? ( "|" + anno.seed_md_tag ) : "" )
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
