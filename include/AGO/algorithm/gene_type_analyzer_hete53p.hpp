#pragma once
#include <cstdlib>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerHete53p
{

  public:

    std::vector< std::map< std::string, std::map< std::size_t, double >>> hete_5p_tables;
    //     smp                miR                 start / end    ppm
    std::vector< std::map< std::string, std::map< std::size_t, double >>> hete_3p_tables;

    GeneTypeAnalyzerHete53p()
    {}

    static void make_hete_table(
            std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            std::vector< std::map< std::string, std::map< std::size_t, double >>>& hete_5p_tables,
            std::vector< std::map< std::string, std::map< std::size_t, double >>>& hete_3p_tables,
            auto& genome_table,
            const std::string& token = ""
            )
    {
        std::size_t arm;
        std::string gene_name;

        bool is_arms = token == "3p" || token == "5p" ? true : false;
        bool is_lens = token != "" && !is_arms ? true : false;

        hete_5p_tables = std::vector< std::map< std::string, std::map< std::size_t, double >>>( bed_samples.size() );
        hete_3p_tables = std::vector< std::map< std::string, std::map< std::size_t, double >>>( bed_samples.size() );

        std::size_t* end_5p;
        std::size_t* end_3p;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.annotation_info_.empty() || raw_bed.annotation_info_[0].empty() ) continue;
                if( biotype == "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != "miRNA" && raw_bed.annotation_info_[0][0] != "mirtron" ) continue;
                if( biotype != "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != biotype ) continue;

                for( int i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    gene_name = raw_bed.annotation_info_[0][ i+1 ];
                    arm = std::stoi( GeneTypeAnalyzerCounting::get_arm( gene_name ).substr( 0, 1 ));

                    if( is_arms && token != ( std::to_string( arm ) + "p"     )) continue;
                    if( is_lens && token !=   std::to_string( raw_bed.length_ )) continue;

                    switch( raw_bed.strand_ )
                    {
                        case '+' :
                            end_5p = &raw_bed.start_;
                            end_3p = &raw_bed.end_;
                            break;

                        case '-' :
                            end_5p = &raw_bed.end_;
                            end_3p = &raw_bed.start_;
                            break;
                    }

                    if( hete_5p_tables[ smp ][ gene_name ].find( *end_5p ) == hete_5p_tables[ smp ][ gene_name ].end() )
                        hete_5p_tables[ smp ][ gene_name ][ *end_5p ] = 0;

                    if( hete_3p_tables[ smp ][ gene_name ].find( *end_3p ) == hete_3p_tables[ smp ][ gene_name ].end() )
                        hete_3p_tables[ smp ][ gene_name ][ *end_3p ] = 0;

                    hete_5p_tables[ smp ][ gene_name ][ *end_5p ] += raw_bed.ppm_;
                    hete_3p_tables[ smp ][ gene_name ][ *end_3p ] += raw_bed.ppm_;
                }
            }

            format_hete_table( hete_5p_tables[ smp ]);
            format_hete_table( hete_3p_tables[ smp ]);
        }
    }

    static void format_hete_table( std::map< std::string, std::map< std::size_t, double >>& hete_tables )
    {
        double ppm_sum;
        std::size_t max_end;

        for( auto& anno : hete_tables )
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

    static void output_heterorgeneity(
            const std::string& output_path,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::map< std::string, std::map< std::size_t, double >>>& hete_tables,
            const std::string& token
            )
    {
        std::vector< std::string > split;
        std::set< std::string > anno_idx;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            if( split.size() > 2 ) for( std::size_t i = 1; i < split.size() -2; ++i )
                    split[0] = split[0] + "_" + split[i];

            anno_idx.emplace( split[0] );
        }

        std::ofstream output( output_path + "heterorgeneity_" + token + ".text" );
        output << "heterorgeneity-" << token;

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : anno_idx )
        {
            output << "\n" << anno;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( hete_tables[ smp ][ anno ].find( 0 ) != hete_tables[ smp ][ anno ].end() )
                {
                    output << "\t" << hete_tables[ smp ][ anno ][0];
                }
                else output << "\t0";
            }
        }

        output.close();
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
