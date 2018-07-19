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
        std::size_t len = 0;
        std::size_t max_len = 0;

        double sum = 0.0;

        isT2C = false;
        idx_string = "";
        md_string = "";
        tc_string = "";

        std::vector< std::map< std::size_t, double >> md_vecs( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> tc_vecs( bed_samples.size(), std::map< std::size_t, double >() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                len = (int)(raw_bed.length_) - (int)(raw_bed.tail_length_);

                if( len > max_len ) max_len = len;

                // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
                {
                    std::size_t i = 0; // do first priority
                    if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                    {
                        if(( raw_bed.annotation_info_[i][0] == biotype ) ||
                           ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))) 
                        {
                            for( int j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                            {
                                if( raw_bed.md_map.size() != 0 )
                                {
                                    for( auto& pos : raw_bed.md_map )
                                    {
                                        if( md_vecs[ smp ].find( pos.first +1 ) == md_vecs[ smp ].end() )
                                            md_vecs[ smp ][ pos.first +1 ] = 0.0;

                                        md_vecs[ smp ][ pos.first +1 ] += raw_bed.ppm_;
                                    }
                                }

                                if( raw_bed.tc_set.size() != 0 )
                                {
                                    isT2C = true;
                                    for( auto& pos : raw_bed.tc_set )
                                    {
                                        if( tc_vecs[ smp ].find( pos +1 ) == tc_vecs[ smp ].end() )
                                            tc_vecs[ smp ][ pos +1 ] = 0.0;

                                        tc_vecs[ smp ][ pos +1 ] += raw_bed.ppm_;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            if( smp != 0 )
            {
                md_string += ",\n";
                tc_string += ",\n";
            }

            md_string += "                        [ 'MD-" + bed_samples[ smp ].first + "'";
            tc_string += "                        [ 'TC-" + bed_samples[ smp ].first + "'";

            for( std::size_t pos = 1; pos <= max_len; ++pos )
            {
                if( md_vecs[ smp ].find( pos ) != md_vecs[ smp ].end() )
                    md_string += ", " + std::to_string( md_vecs[ smp ][ pos ] );
                else
                    md_string += ", 0";

                if( tc_vecs[ smp ].find( pos ) != tc_vecs[ smp ].end() )
                    tc_string += ", " + std::to_string( tc_vecs[ smp ][ pos ] );
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
        std::ofstream output( output_name + "index.html" );

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
