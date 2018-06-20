#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerRanking
{
  public:

    GeneTypeAnalyzerRanking()
    {}

    static void output_ranking(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            const std::string& token
            )
    {
        std::ofstream output( output_name + token + ".tsv" );
        output << "Annotation";

        std::string anno_name;
        std::size_t index = 0;

        std::vector< std::string > anno_names;
        std::vector< std::vector< std::pair< std::size_t, double >>> rank_temps( anno_table_tail.size() );

        double value;
        double gm = 0.0;
        double pm = 0.0;

        for( auto& anno : ano_len_idx.first )
        {
            anno_name = anno + ( anno_mark[0].find( anno ) != anno_mark[0].end() ? anno_mark[0][ anno ] : "" );
            anno_names.emplace_back( anno_name );

            for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
            {
                gm = 0.0;
                pm = 0.0;

                if( token != "PM" )
                {
                    if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                                gm += anno_table_tail[ smp ][5][ anno ][ len ];
                        }
                }

                if( token != "GM" )
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

                value = ( token == "GMPM" ? gm + pm : ( token == "GM" ? gm : ( token == "PM" ? pm : (( gm + pm ) < 1 ? 0 : ( pm * 100 / ( gm + pm ))))));
                rank_temps[ smp ].emplace_back( std::make_pair( index, value ));
            }

            index++;
        }

        sort_pair( rank_temps );
        make_rank( rank_temps );
        reordring( rank_temps );

        for( auto& smp : bed_samples ) output << "\t" << smp.first;
        for( std::size_t idx = 0; idx < rank_temps[0].size(); ++idx )
        {
            output << "\n" << anno_names[ rank_temps[0][ idx ].first ] << std::setprecision( 0 ) << std::fixed;
            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                output << "\t" << rank_temps[ smp ][ idx ].second;
        }

        output << "\n";
        output.close();
    }

    static void sort_pair( std::vector< std::vector< std::pair< std::size_t, double >>>& pair_vecs )
    {
        for( auto& pairs : pair_vecs )
        {
            std::sort( pairs.begin(), pairs.end(), [] ( const std::pair< std::size_t, double >& a, const std::pair< std::size_t, double >& b )
            { 
                if( a.second == b.second )
                    return a.first < b.first;
                else
                    return a.second < b.second;
            });
        }
    }

    static void make_rank( std::vector< std::vector< std::pair< std::size_t, double >>>& pair_vecs )
    {
        double rank;
        for( auto& pairs : pair_vecs )
        {
            rank = 0;
            for( auto& pair : pairs )
            {
                rank++;
                pair.second = rank;
            }
        }
    }

    static void reordring( std::vector< std::vector< std::pair< std::size_t, double >>>& pair_vecs )
    {
        double sum;
        std::map< std::size_t, double > temp;
        std::vector< std::pair< std::size_t, double >> sum_rank;
        std::vector< std::vector< std::pair< std::size_t, double >>> result( pair_vecs.size() );

        for( std::size_t smp = 0; smp < pair_vecs.size(); ++smp )
        {
            temp.clear();
            for( auto& pair : pair_vecs[ smp ]) temp.emplace( pair );

            pair_vecs[ smp ].clear();
            for( auto& pair : temp ) pair_vecs[ smp ].emplace_back( pair );
        }

        for( std::size_t anno = 0; anno < pair_vecs[ 0 ].size(); ++anno )
        {
            sum = 0;
            for( std::size_t samp = 0; samp < pair_vecs.size(); ++samp )
                sum += pair_vecs[ samp ][ anno ].second;
            sum_rank.emplace_back( std::make_pair( anno, sum ));
        }

        std::sort( sum_rank.begin(), sum_rank.end(), [] ( const std::pair< std::size_t, double >& a, const std::pair< std::size_t, double >& b )
        { 
            if( a.second == b.second )
                return a.first > b.first;
            else
                return a.second > b.second;
        });

        for( auto& rank : sum_rank )
        {
            for( std::size_t samp = 0; samp < pair_vecs.size(); ++samp )
                result[ samp ].emplace_back( pair_vecs[ samp ][ rank.first ] );
        }

        pair_vecs = result;
    }
};

} // end of namespace algorithm
} // end of namespace ago
