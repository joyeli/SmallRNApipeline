#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDifference
{
  public:

    GeneTypeAnalyzerDifference()
    {}

    double get_sum( const std::vector< double >& vec )
    {
        double res = 0;
        for( auto& val : vec )
            res += val;
        return res;
    }

    std::tuple< double, std::size_t, std::size_t > get_difference( const std::vector< double >& vec )
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

            if( max.first < vec[ smp ] )
            {
                max = std::make_pair( vec[ smp ], smp );
            }

            if( min.first > vec[ smp ] )
            {
                min = std::make_pair( vec[ smp ], smp );
            }
        }

        std::vector< double > folds;
        double fold = (( max.first/sum )-( min.first/sum ))/( min.first/sum );

        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == max.second || smp == min.second ) continue;
            folds.emplace_back( (( vec[ smp ]/sum )-( min.first/sum ))/( min.first/sum ));
        }

        bool fold_tag = false;

        for( auto& fd : folds )
            if( fd >= 1 ) fold_tag = true;

        return{ fold_tag ? ( -1 * fold ) : fold, max.second, min.second };
    }

    std::vector< std::size_t > get_lengths( const std::vector< std::map< std::size_t, double >>& len_vec )
    {
        double max;
        std::size_t len;
        std::vector< std::size_t > res( len_vec.size(), std::size_t() );

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            len = 0;
            max = 0.0;

            for( auto& lens : len_vec[ smp ] )
            {
                if( max < lens.second )
                {
                    len = lens.first;
                    max = lens.second;
                }
            }

            res[ smp ] = len;
        }

        return res;
    }

    std::tuple< std::size_t, std::size_t, std::size_t > get_length_difference( const std::vector< std::size_t >& len_vec )
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

            if( max.first < len_vec[ smp ] )
            {
                max = std::make_pair( len_vec[ smp ], smp );
            }

            if( min.first > len_vec[ smp ] )
            {
                min = std::make_pair( len_vec[ smp ], smp );
            }
        }

        return { max.first - min.first, max.second, min.second };
    }

    std::vector< std::pair< std::string, double >> get_arms( const std::vector< std::map< std::string, double >>& arm_vec )
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

    std::string get_arm_difference( const std::vector< std::pair< std::string, double >>& arm_vec )
    {
        std::string res = "Y";
        std::set< std::string > samp;
        for( auto& arm : arm_vec ) samp.emplace( arm.first );
        if( samp.size() == 1 ) res = "N";
        return res;
    }

    void sort_difference( std::vector< std::pair< double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
            { return a.first > b.first; });
    }

    void output_loading_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables,
            const double& sudo_count
            )
    {
        std::ofstream out_loading( output_path + "/LoadingDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/LoadingDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";
        out_tailing << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";

        std::vector< double > vec_loading;
        std::vector< double > vec_tailing;

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::tuple< double, std::size_t, std::size_t > df_tuple_loading;
        std::tuple< double, std::size_t, std::size_t > df_tuple_tailing;

        double pm = 0.0;
        double gmpm = 0.0;
        double sum_loading = 0.0;

        std::string res_loading;
        std::string res_tailing;

        for( auto& anno : ano_len_idx.first )
        {
            vec_loading.clear();
            vec_tailing.clear();

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                pm   = counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];
                gmpm = counting_tables[ smp ].first[ anno ][ 0 ];

                vec_loading.emplace_back( gmpm );
                vec_tailing.emplace_back( pm < 1 ? sudo_count : ( pm / gmpm * 100 ));
            }

            sum_loading = get_sum( vec_loading );

            df_tuple_loading = get_difference( vec_loading );
            df_tuple_tailing = get_difference( vec_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_loading += std::to_string( std::get<0>( df_tuple_loading )) + "\t" + bed_samples[ std::get<1>( df_tuple_loading )].first + ":" + bed_samples[ std::get<2>( df_tuple_loading )].first;
            res_tailing += std::to_string( std::get<0>( df_tuple_tailing )) + "\t" + bed_samples[ std::get<1>( df_tuple_tailing )].first + ":" + bed_samples[ std::get<2>( df_tuple_tailing )].first;

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }

    void output_length_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_loading( output_path + "/LengthDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/LengthDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";
        out_tailing << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::vector< std::size_t > len_loading;
        std::vector< std::size_t > len_tailing;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::size_t, std::tuple< double, double >>>> smp_anno_len_gmpms;

        std::vector< std::map< std::size_t, double >> vec_len_loading( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_tailing( bed_samples.size(), std::map< std::size_t, double >() );

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_loading;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_tailing;

        std::string res_loading;
        std::string res_tailing;

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

                /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] ) += counting_tables[ smp ].first[ anno ][ 0 ];

                for( auto& len : ano_len_idx.second )
                {
                    if( smp_anno_len_gmpms[ smp ][ anno_tag ].find( len ) == smp_anno_len_gmpms[ smp ][ anno_tag ].end() )
                        smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] = { 0, 0 };

                    /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) +=   counting_tables[ smp ].first[ anno ][ len ];
                    /* PM */ std::get<1>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) += ( counting_tables[ smp ].first[ anno ][ len ] * counting_tables[ smp ].second[ anno ][ len ] );
                }
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_len_loading = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );
            vec_len_tailing = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& len : smp_anno_len_gmpms[ smp ][ anno ] )
                {
                    if( len.first == 0 )
                    {
                        vec_loading.emplace_back( std::get<0>( len.second ) );
                        continue;
                    }

                    vec_len_loading[ smp ][ len.first ] = std::get<0>( len.second );
                    vec_len_tailing[ smp ][ len.first ] = std::get<1>( len.second );
                }
            }

            sum_loading = get_sum( vec_loading );

            len_loading = get_lengths( vec_len_loading );
            len_tailing = get_lengths( vec_len_tailing );

            df_tuple_loading = get_length_difference( len_loading );
            df_tuple_tailing = get_length_difference( len_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_loading += std::to_string( std::get<0>( df_tuple_loading )) + "\t";
            res_tailing += std::to_string( std::get<0>( df_tuple_tailing )) + "\t";

            res_loading += bed_samples[ std::get<1>( df_tuple_loading )].first + ":" + bed_samples[ std::get<2>( df_tuple_loading )].first + "\t";
            res_tailing += bed_samples[ std::get<1>( df_tuple_tailing )].first + ":" + bed_samples[ std::get<2>( df_tuple_tailing )].first + "\t";

            res_loading += std::to_string( len_loading[ std::get<1>( df_tuple_loading )]) + ":" + std::to_string( len_loading[ std::get<2>( df_tuple_loading )]);
            res_tailing += std::to_string( len_tailing[ std::get<1>( df_tuple_tailing )]) + ":" + std::to_string( len_tailing[ std::get<2>( df_tuple_tailing )]);

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }

    void output_arms_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_loading( output_path + "/ArmDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/ArmDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tisLoadingDifference";
        out_tailing << "Annotation\tTotal\tisLoadingDifference";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_loading << "\tSample1";
            else out_loading << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_tailing << "\tArm1";
            else out_tailing << ":Arm" << smp + 1 ;

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::string, std::tuple< double, double >>>> smp_anno_arm_gmpms;

        std::vector< std::map< std::string, double >> vec_arm_loading( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::map< std::string, double >> vec_arm_tailing( bed_samples.size(), std::map< std::string, double >() );

        std::vector< std::pair< std::string, double >> arm_loading;
        std::vector< std::pair< std::string, double >> arm_tailing;

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::string is_diff_loading;
        std::string is_diff_tailing;

        std::string res_loading;
        std::string res_tailing;

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

                /*GMPM*/ std::get<0>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += counting_tables[ smp ].first[ anno ][ 0 ];
                /* PM */ std::get<1>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_arm_loading = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );
            vec_arm_tailing = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& arm : smp_anno_arm_gmpms[ smp ][ anno ] )
                {
                    vec_loading.emplace_back( std::get<0>( arm.second ) );
                    vec_arm_loading[ smp ][ arm.first ] =   std::get<0>( arm.second );
                    vec_arm_tailing[ smp ][ arm.first ] = ( std::get<1>( arm.second ) < 1 ? 0 : std::get<1>( arm.second )) / std::get<0>( arm.second ) * 100;
                }
            }

            sum_loading = get_sum( vec_loading );

            arm_loading = get_arms( vec_arm_loading );
            arm_tailing = get_arms( vec_arm_tailing );

            is_diff_loading = get_arm_difference( arm_loading );
            is_diff_tailing = get_arm_difference( arm_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_loading;
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_tailing;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_loading += "\t" + bed_samples[ smp ].first;
                else res_loading += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_loading += "\t" + arm_loading[ smp ].first;
                else res_loading += ":" + arm_loading[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tailing += "\t" + bed_samples[ smp ].first;
                else res_tailing += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tailing += "\t" + arm_tailing[ smp ].first;
                else res_tailing += ":" + arm_tailing[ smp ].first;

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
