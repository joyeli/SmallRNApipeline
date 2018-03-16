#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDifference
{
  public:

    GeneTypeAnalyzerDifference()
    {}

    static double get_sum( const std::vector< double >& vec )
    {
        double res = 0;
        for( auto& val : vec )
            res += val;
        return res;
    }

    static void sort_difference( std::vector< std::pair< double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
            { return a.first > b.first; });
    }



    static std::tuple< double, std::size_t, std::size_t > get_difference( const std::vector< double >& vec )
    {
        double sum = get_sum( vec );
        std::pair< double, std::size_t > max;
        std::pair< double, std::size_t > min;

        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == 0 )
            {
                max = std::make_pair( vec[ smp ], smp );
                min = std::make_pair( vec[ smp ], smp );
                continue;
            }

            if( max.first < vec[ smp ] ) max = std::make_pair( vec[ smp ], smp );
            if( min.first > vec[ smp ] ) min = std::make_pair( vec[ smp ], smp );
        }

        double max_tmp = ( max.first == 0.0 ? 0.000001 : max.first ) / sum;
        double min_tmp = ( min.first == 0.0 ? 0.000001 : min.first ) / sum;
        double temp = max_tmp - min_tmp;

        double fold = ( sum == 0.0 ? 0.0 : ( temp / min_tmp ));
        std::vector< double > folds;

        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == max.second || smp == min.second ) continue;
            folds.emplace_back( sum == 0 ? 0 : ((( vec[ smp ]/sum )-( min.first/sum ))/( min.first/sum )));
        }

        bool fold_tag = false;
        for( auto& fd : folds ) if( fd >= 1 ) fold_tag = true;
        return{ fold_tag ? ( -1 * fold ) : fold, max.second, min.second };
    }

    static void output_loading_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& token
            )
    {
        std::ofstream output( output_path + "LoadingDifference_" + token + ".tsv" );
        output << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";

        std::vector< double > vec_gmpm;
        std::vector< double > vec_temp;
        std::vector< std::pair< double, std::string >> out_temp;
        std::tuple< double, std::size_t, std::size_t > df_tuple;

        double gm = 0.0;
        double pm = 0.0;
        double sum_loading = 0.0;

        for( auto& anno : ano_len_idx.first )
        {
            vec_gmpm.clear();
            vec_temp.clear();

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                gm = 0.0;
                pm = 0.0;

                if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                    for( auto& len : ano_len_idx.second )
                    {
                        if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                            gm += anno_table_tail[ smp ][5][ anno ][ len ];
                    }

                for( std::size_t i = 0; i < 5; i++ )
                {
                    if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                pm += anno_table_tail[ smp ][i][ anno ][ len ];
                        }
                }

                vec_gmpm.emplace_back( gm + pm );
                vec_temp.emplace_back(
                        ( token == "GMPM" ? ( gm + pm )
                        : token == "GM"   ? gm
                        : token == "PM"   ? pm
                        : ( pm < 1 ? 0.0
                        : ( pm * 100 / ( gm + pm )) ))
                        );
            }

            sum_loading = get_sum( vec_gmpm );
            df_tuple = get_difference( vec_temp );

            out_temp.emplace_back( std::make_pair(
                        sum_loading,
                        "\n" + anno + "\t"
                            + std::to_string( sum_loading ) + "\t"
                            + std::to_string( std::get<0>( df_tuple )) + "\t"
                            + bed_samples[ std::get<1>( df_tuple )].first + ":" + bed_samples[ std::get<2>( df_tuple )].first
                        ));
        }

        sort_difference( out_temp );
        for( auto& temp : out_temp )
            output << temp.second;

        output << "\n";
        output.close();
    }



    static std::vector< std::size_t > get_lengths( const std::vector< std::map< std::size_t, double >>& len_vec )
    {
        double max;
        std::size_t len;
        std::vector< std::size_t > res( len_vec.size(), std::size_t() );

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            len = 0;
            max = 0.0;

            for( auto& lens : len_vec[ smp ] ) if( max < lens.second )
            {
                len = lens.first;
                max = lens.second;
            }

            res[ smp ] = len;
        }

        return res;
    }

    static std::tuple< std::size_t, std::size_t, std::size_t > get_length_difference( const std::vector< std::size_t >& len_vec )
    {
        std::pair< std::size_t, std::size_t > max;
        std::pair< std::size_t, std::size_t > min;

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            if( smp == 0 )
            {
                max = std::make_pair( len_vec[ smp ], smp );
                min = std::make_pair( len_vec[ smp ], smp );
                continue;
            }

            if( max.first < len_vec[ smp ] ) max = std::make_pair( len_vec[ smp ], smp );
            if( min.first > len_vec[ smp ] ) min = std::make_pair( len_vec[ smp ], smp );
        }

        return { max.first - min.first, max.second, min.second };
    }

    static void output_length_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& token
            )
    {
        std::ofstream output( output_path + "LengthDifference_" + token + ".tsv" );
        output << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";

        double sum_loading = 0.0;
        std::vector< double > vec_loading;
        std::vector< std::size_t > lens;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::size_t, std::tuple< double, double >>>> smp_anno_len_gmpms;

        std::vector< std::map< std::size_t, double >> vec_len( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::pair< double, std::string >> out_temp;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple;

        std::vector< std::string > anno_splits; 
        std::string anno_tag;

        for( auto& anno : ano_len_idx.first )
        {
            anno_splits.clear();

            boost::iter_split( anno_splits, anno, boost::algorithm::first_finder( "_" ));
            anno_tag = anno_splits[0];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "_" + anno_splits[i] );
            anno_index.emplace( anno_tag );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( smp_anno_len_gmpms[ smp ][ anno_tag ].find( 0 ) == smp_anno_len_gmpms[ smp ][ anno_tag ].end() )
                    smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] = { 0, 0 };

                for( auto& len : ano_len_idx.second )
                    if( smp_anno_len_gmpms[ smp ][ anno_tag ].find( len ) == smp_anno_len_gmpms[ smp ][ anno_tag ].end() )
                        smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] = { 0, 0 };

                if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                    for( auto& len : ano_len_idx.second )
                    {
                        if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                        {
                            /*GM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) += anno_table_tail[ smp ][5][ anno ][ len ];
                            /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] ) += anno_table_tail[ smp ][5][ anno ][ len ];
                        }
                    }

                for( std::size_t i = 0; i < 5; i++ )
                {
                    if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                            {
                                /*PM*/ std::get<1>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) += anno_table_tail[ smp ][i][ anno ][ len ];
                                /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] ) += anno_table_tail[ smp ][i][ anno ][ len ];
                            }
                        }
                }
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_len = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& len : smp_anno_len_gmpms[ smp ][ anno ] )
                {
                    if( len.first == 0 )
                    {
                        vec_loading.emplace_back( std::get<0>( len.second ) );
                        continue;
                    }

                    vec_len[ smp ][ len.first ] =
                             ( token == "GMPM" ? ( std::get<0>( len.second ) + std::get<1>( len.second ))
                             : token == "GM"   ? ( std::get<0>( len.second ) )
                             : token == "PM"   ? ( std::get<1>( len.second ) )
                             : ( std::get<1>( len.second ) < 1 ? 0.0
                             : ( std::get<1>( len.second ) * 100 / ( std::get<0>( len.second ) + std::get<1>( len.second )) )) 
                             );
                }
            }

            sum_loading = get_sum( vec_loading );

            lens = get_lengths( vec_len );
            df_tuple = get_length_difference( lens );

            out_temp.emplace_back( std::make_pair(
                        sum_loading,
                        "\n" + anno + "\t"
                            + std::to_string( sum_loading ) + "\t"
                            + std::to_string( std::get<0>( df_tuple )) + "\t"
                            + bed_samples[ std::get<1>( df_tuple )].first + ":" + bed_samples[ std::get<2>( df_tuple )].first + "\t"
                            + std::to_string( lens[ std::get<1>( df_tuple )]) + ":" + std::to_string( lens[ std::get<2>( df_tuple )])
                        ));
        }

        sort_difference( out_temp );
        for( auto& temp : out_temp )
            output << temp.second;

        output << "\n";
        output.close();
    }



    static std::vector< std::pair< std::string, double >> get_arms( const std::vector< std::map< std::string, double >>& arm_vec )
    {
        std::vector< std::pair< std::string, double >> res;
        std::pair< std::string, double > tmp = std::make_pair( "", 0.0 );

        for( std::size_t smp = 0; smp < arm_vec.size(); ++smp )
        {
            for( auto& arm : arm_vec[ smp ] )
            {
                if( tmp.second == 0.0 ) tmp = arm;
                if( tmp.second < arm.second ) tmp = arm;
            }

            res.emplace_back( tmp );
            tmp = std::make_pair( "", 0.0 );
        }

        return res;
    }

    static std::string get_arm_difference( const std::vector< std::pair< std::string, double >>& arm_vec )
    {
        std::string res = "Y";
        std::set< std::string > samp;
        for( auto& arm : arm_vec ) samp.emplace( arm.first );
        if( samp.size() == 1 ) res = "N";
        return res;
    }

    static void output_arms_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& token
            )
    {
        std::ofstream output( output_path + "ArmDifference_" + token + ".tsv" );
        output << "Annotation\tTotal\tisLoadingDifference";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) output << "\tSample1";
            else output << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) output << "\tArm1";
            else output << ":Arm" << smp + 1 ;

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::string, std::tuple< double, double >>>> smp_anno_arm_gmpms;

        std::vector< std::map< std::string, double >> vec_arm( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::pair< std::string, double >> arms;
        std::vector< std::pair< double, std::string >> out_temp;

        std::string is_diff;
        std::string res;

        std::vector< std::string > anno_splits; 
        std::string anno_tag;
        std::string anno_arm;

        for( auto& anno : ano_len_idx.first )
        {
            anno_splits.clear();

            boost::iter_split( anno_splits, anno, boost::algorithm::first_finder( "_" ));
            anno_tag = anno_splits[0];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "_" + anno_splits[i] );
            anno_splits.clear();

            boost::iter_split( anno_splits, anno_tag, boost::algorithm::first_finder( "-" ));

            anno_tag = anno_splits[0];
            anno_arm = anno_splits[ anno_splits.size() -1 ];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "-" + anno_splits[i] );
            anno_index.emplace( anno_tag );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( smp_anno_arm_gmpms[ smp ][ anno_tag ].find( anno_arm ) == smp_anno_arm_gmpms[ smp ][ anno_tag ].end() )
                    smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] = { 0, 0 };

                if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                    for( auto& len : ano_len_idx.second )
                    {
                        if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                            /*GM*/ std::get<0>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += anno_table_tail[ smp ][5][ anno ][ len ];
                    }

                for( std::size_t i = 0; i < 5; i++ )
                {
                    if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                /*PM*/ std::get<1>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += anno_table_tail[ smp ][i][ anno ][ len ];
                        }
                }
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_arm = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& arm : smp_anno_arm_gmpms[ smp ][ anno ] )
                {
                    vec_loading.emplace_back( std::get<0>( arm.second ) + std::get<1>( arm.second ));

                    vec_arm[ smp ][ arm.first ] =
                             ( token == "GMPM" ? ( std::get<0>( arm.second ) + std::get<1>( arm.second ))
                             : token == "GM"   ? ( std::get<0>( arm.second ) )
                             : token == "PM"   ? ( std::get<1>( arm.second ) )
                             : ( std::get<1>( arm.second ) < 1 ? 0.0
                             : ( std::get<1>( arm.second ) * 100 / ( std::get<0>( arm.second ) + std::get<1>( arm.second )) )) 
                             );
                }
            }

            sum_loading = get_sum( vec_loading );

            arms = get_arms( vec_arm );
            is_diff = get_arm_difference( arms );

            res = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res += "\t" + bed_samples[ smp ].first;
                else res += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res += "\t" + arms[ smp ].first;
                else res += ":" + arms[ smp ].first;

            out_temp.emplace_back( std::make_pair( sum_loading, res ));
        }

        sort_difference( out_temp );
        for( auto& temp : out_temp )
            output << temp.second;

        output << "\n";
        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
