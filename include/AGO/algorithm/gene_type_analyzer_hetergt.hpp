#pragma once
#include <cstdlib>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerHetergt
{

  public:

    GeneTypeAnalyzerHetergt()
    {}

    static std::vector< std::map< std::string, std::tuple< double, double, double >>> make_hete_table(
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::string& biotype,
            const std::string& token = ""
            )
    {
        bool is_lens = token != "" ? true : false;
        std::vector< std::map< std::string, std::tuple< double, double, double >>> results( bed_samples.size() );
        //                                              5p      3p      tail

        std::vector< std::vector< 
            std::map< std::string, std::map< std::size_t, double >> // 5p / 3p / tail
        >> hete_tables( bed_samples.size(), std::vector< std::map< std::string, std::map< std::size_t, double >>>( 3 ));

        std::string gene_name;

        std::size_t* end_5p;
        std::size_t* end_3p;
        std::size_t  tailed;

        for( auto& anno : ano_len_idx.first )
        {
            if( results[0].find( anno ) == results[0].end() )
            {
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                    results[ smp ][ anno ] = { 0.0, 0.0, 0.0 };
            }
        }

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                if( raw_bed.annotation_info_.empty() || raw_bed.annotation_info_[0].empty() ) continue;
                if( biotype == "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != "miRNA" && raw_bed.annotation_info_[0][0] != "mirtron" ) continue;
                if( biotype != "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != biotype ) continue;

                for( int i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    if( is_lens && token != std::to_string( raw_bed.length_ - raw_bed.tail_length_ ))
                        continue;

                    gene_name = raw_bed.annotation_info_[0][ i+1 ];

                    switch( raw_bed.strand_ )
                    {
                        case '+' :
                            end_5p = &raw_bed.start_;
                            end_3p = &raw_bed.end_;
                            tailed = *end_3p + raw_bed.tail_length_;
                            break;

                        case '-' :
                            end_5p = &raw_bed.end_;
                            end_3p = &raw_bed.start_;
                            tailed = *end_3p - raw_bed.tail_length_;
                            break;
                    }

                    if( hete_tables[ smp ][0][ gene_name ].find( *end_5p ) == hete_tables[ smp ][0][ gene_name ].end() )
                        hete_tables[ smp ][0][ gene_name ][ *end_5p ] = 0;

                    if( hete_tables[ smp ][1][ gene_name ].find( *end_3p ) == hete_tables[ smp ][1][ gene_name ].end() )
                        hete_tables[ smp ][1][ gene_name ][ *end_3p ] = 0;

                    if( hete_tables[ smp ][2][ gene_name ].find(  tailed ) == hete_tables[ smp ][2][ gene_name ].end() )
                        hete_tables[ smp ][2][ gene_name ][  tailed ] = 0;

                    hete_tables[ smp ][0][ gene_name ][ *end_5p ] += raw_bed.ppm_;
                    hete_tables[ smp ][1][ gene_name ][ *end_3p ] += raw_bed.ppm_;
                    hete_tables[ smp ][2][ gene_name ][  tailed ] += raw_bed.ppm_;
                }
            }

            for( std::size_t hete = 0; hete < hete_tables[ smp ].size(); ++hete )
            {
                format_hete_table( hete_tables[ smp ][ hete ] );

                for( auto& anno : hete_tables[ smp ][ hete ] ) switch( hete )
                {
                    case 0 : std::get<0>( results[ smp ][ anno.first ] ) = anno.second[0]; break;
                    case 1 : std::get<1>( results[ smp ][ anno.first ] ) = anno.second[0]; break;
                    case 2 : std::get<2>( results[ smp ][ anno.first ] ) = anno.second[0]; break;
                }
            }
        }

        return results;
    }

    static void format_hete_table( std::map< std::string, std::map< std::size_t, double >>& hete_table )
    {
        double ppm_sum;
        std::size_t max_end;

        for( auto& anno : hete_table )
        {
            max_end = get_max_end( anno.second );
            ppm_sum = get_ppm_sum( anno.second );

            for( auto& end : anno.second )
                end.second = (double)(std::abs( (int)(max_end) - (int)(end.first) )) * end.second;

            anno.second[0] = get_ppm_sum( anno.second ) / ppm_sum;
        }
    }

    static std::size_t get_max_end( std::map< std::size_t, double >& ends_table )
    {
        std::size_t max_end;
        double max_ppm = 0.0;

        for( auto& end : ends_table )
        {
            if( end.second > max_ppm )
            {
                max_end = end.first;
                max_ppm = end.second;
            }
        }

        return max_end;
    }

    static double get_ppm_sum( std::map< std::size_t, double >& ends_table )
    {
        double ppm_sum = 0.0;
        for( auto& end : ends_table ) ppm_sum += end.second;
        return ppm_sum;
    }

    static void debug( std::map< std::string, std::map< std::size_t, double >>& hete_table )
    {
        std::ofstream output( "./debug.tsv" );

        for( auto& gene : hete_table )
            for( auto& end : gene.second )
                output << gene.first << "\t" << end.first << "\t" << end.second << "\n";
        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
