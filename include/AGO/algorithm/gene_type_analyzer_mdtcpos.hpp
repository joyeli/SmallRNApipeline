#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerMDTCpos
{
    static bool isT2C;
    static std::string idx_string;
    static std::string md_string;
    static std::string tc_string;

  public:

    GeneTypeAnalyzerMDTCpos()
    {}

    static void make_mdtcpos( std::vector< BedSampleType >& bed_samples, const std::string& biotype )
    {
        double sum = 0;
        std::size_t len = 0;
        std::size_t max_len = 0;

        isT2C = false;
        idx_string = "";
        md_string = "";
        tc_string = "";

        std::vector< std::map< std::size_t, double >> md_vecs( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> tc_vecs( bed_samples.size(), std::map< std::size_t, double >() );

        for( std::size_t idx = 0; idx < bed_samples.size(); ++idx )
        {
            sum = 0;

            for( auto& raw_bed : bed_samples[ idx ].second )
            {
                for( auto& raw_bed_info : raw_bed.annotation_info_ )
                {
                    for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                    {
                        if( raw_bed_info[i] != biotype ) continue;

                        sum += ( raw_bed.reads_count_ / raw_bed.multiple_alignment_site_count_ );
                        len = ( int )raw_bed.length_ - ( int )raw_bed.tail_length_;

                        if( len > max_len ) max_len = len;

                        if( raw_bed.md_map.size() != 0 )
                        {
                            for( auto& pos : raw_bed.md_map )
                            {
                                if( md_vecs[ idx ].find( pos.first +1 ) == md_vecs[ idx ].end() )
                                    md_vecs[ idx ][ pos.first +1 ] = 0;

                                md_vecs[ idx ][ pos.first +1 ]++;
                            }
                        }

                        if( raw_bed.tc_set.size() != 0 )
                        {
                            isT2C = true;
                            for( auto& pos : raw_bed.tc_set )
                            {
                                if( tc_vecs[ idx ].find( pos +1 ) == tc_vecs[ idx ].end() )
                                    tc_vecs[ idx ][ pos +1 ] = 0;

                                tc_vecs[ idx ][ pos +1 ]++;
                            }
                        }
                    }
                }
            }

            for( auto& pos : md_vecs[ idx ] ) pos.second = pos.second / sum;
            for( auto& pos : tc_vecs[ idx ] ) pos.second = pos.second / sum;
        }

        for( std::size_t idx = 0; idx < bed_samples.size(); ++idx )
        {
            if( idx != 0 )
            {
                md_string += ",\n";
                tc_string += ",\n";
            }

            md_string += "                        [ 'MD-" + bed_samples[ idx ].first + "'";
            tc_string += "                        [ 'TC-" + bed_samples[ idx ].first + "'";

            for( std::size_t pos = 1; pos <= max_len; ++pos )
            {
                if( md_vecs[ idx ].find( pos ) != md_vecs[ idx ].end() )
                    md_string += ", " + std::to_string( md_vecs[ idx ][ pos ] );
                else
                    md_string += ", 0";

                if( tc_vecs[ idx ].find( pos ) != tc_vecs[ idx ].end() )
                    tc_string += ", " + std::to_string( tc_vecs[ idx ][ pos ] );
                else
                    tc_string += ", 0";
            }

            md_string += " ]";
            tc_string += " ]";
        }

        for( std::size_t pos = 1; pos <= max_len; ++pos )
            idx_string += "'" + std::to_string( pos ) + "'" + ( pos != max_len ? ", " : "" );
    }

    static void output_mdtcpos_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "MDTCpos.html" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <head>" << "\n";
        output << "        <link href='https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.2/c3.min.css' rel='stylesheet'>" << "\n";
        output << "        <script src='https://cdnjs.cloudflare.com/ajax/libs/d3/4.13.0/d3.js'></script>" << "\n";
        output << "        <script src='https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.2/c3.min.js'></script>" << "\n";
        output << "    </head>" << "\n";
        output << "    <body>" << "\n";
        output << "        <div id='chart'></div>" << "\n";
        output << "        <script>" << "\n";
        output << "            var chart =  c3.generate({" << "\n";
        output << "                bindto: '#chart'," << "\n";
        output << "                size: {" << "\n";
        output << "                    height: 780" << "\n";
        output << "                }," << "\n";
        output << "                data: {" << "\n";
        output << "                    columns: [" << "\n";
        output << md_string + ( !isT2C ? "" : ( ",\n" + tc_string )) << "\n";
        output << "                    ]," << "\n";
        output << "                    type: 'spline'" << "\n";
        output << "                }," << "\n";
        output << "                axis: {" << "\n";
        output << "                    x: {" << "\n";
        output << "                        type: 'category'," << "\n";
        output << "                        categories: [" << "\n";
        output << "                            " << idx_string << "\n";
        output << "                        ]" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "            });" << "\n";
        output << "        </script>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

bool GeneTypeAnalyzerMDTCpos::isT2C = false;
std::string GeneTypeAnalyzerMDTCpos::idx_string = "";
std::string GeneTypeAnalyzerMDTCpos::md_string = "";
std::string GeneTypeAnalyzerMDTCpos::tc_string = "";

} // end of namespace algorithm
} // end of namespace ago
