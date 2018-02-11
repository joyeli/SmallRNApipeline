#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerValplot
{

  public:

    GeneTypeAnalyzerValplot()
    {}

    static void output_valplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >> anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            const std::string& token
            )
    {
        std::ofstream output( output_name + token + ".tsv" );
        output << "Annotation";

        double gm = 0.0;
        double pm = 0.0;

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : ano_len_idx.first )
        {
            output << "\n" << anno << std::setprecision( 0 ) << std::fixed
                << ( anno_mark[0].find( anno ) != anno_mark[0].end() ? anno_mark[0][ anno ] : "" );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                gm = 0.0;
                pm = 0.0;

                if( token != "PM" || token == "Tailing" )
                {
                    if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                                gm += anno_table_tail[ smp ][5][ anno ][ len ];
                        }
                }

                if( token != "GM" || token == "Tailing" )
                {
                    for( std::size_t i = 0; i < 5; i++ )
                    {
                        if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                            for( auto& len : ano_len_idx.second )
                            {
                                if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                    pm += anno_table_tail[ smp ][i][ anno ][ len ];
                            }
                    }
                }

                output << "\t" << ( token == "GMPM" ? gm + pm : ( token == "GM" ? gm : ( token == "PM" ? pm : ( pm * 100 / ( gm + pm )))));
            }
        }

        output << "\n";
        output.close();
    }

    static void output_valplot_visualization( const std::string& output_name )
    {
        // std::ofstream output( output_name + "index.php" );
        // output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
