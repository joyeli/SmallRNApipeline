#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDifference
{
    std::string output_path;
    std::string difference;

  public:

    GeneTypeAnalyzerDifference()
        : output_path( "" )
        , difference( "/Difference/" )
    {}

    GeneTypeAnalyzerDifference(
            const std::string& biotype,
            const std::string& output_path_,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
        : output_path( output_path_ + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" ))
        , difference( "/Difference/" )
    {
        boost::filesystem::create_directory( boost::filesystem::path( output_path + biotype + difference ));

        output_loading_difference( output_path + biotype + difference, bed_samples, ano_len_idx, anno_table_tail );
        output_length_difference(  output_path + biotype + difference, bed_samples, ano_len_idx, anno_table_tail );
        if( biotype == "miRNA" ) output_arms_difference( output_path + biotype + difference, bed_samples, ano_len_idx, anno_table_tail );
    }

    double get_sum( const std::vector< double >& vec )
    {
        double res = 0;
        for( auto& val : vec )
            res += val;
        return res;
    }

    void sort_difference( std::vector< std::pair< double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
            { return a.first > b.first; });
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

            if( max.first < vec[ smp ] ) max = std::make_pair( vec[ smp ], smp );
            if( min.first > vec[ smp ] ) min = std::make_pair( vec[ smp ], smp );
        }

        std::vector< double > folds;
        double fold = ( sum == 0 ? 0 : ((( max.first/sum )-( min.first/sum ))/( min.first/sum )));

        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == max.second || smp == min.second ) continue;
            folds.emplace_back( sum == 0 ? 0 : ((( vec[ smp ]/sum )-( min.first/sum ))/( min.first/sum )));
        }

        bool fold_tag = false;
        for( auto& fd : folds ) if( fd >= 1 ) fold_tag = true;
        return{ fold_tag ? ( -1 * fold ) : fold, max.second, min.second };
    }

    void output_loading_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        std::ofstream out_gmpm( output_path + "LoadingDifference_GMPM.tsv"    );
        std::ofstream out_gm  ( output_path + "LoadingDifference_GM.tsv"      );
        std::ofstream out_pm  ( output_path + "LoadingDifference_PM.tsv"      );
        std::ofstream out_tail( output_path + "LoadingDifference_Tailing.tsv" );

        out_gmpm << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";
        out_gm   << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";
        out_pm   << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";
        out_tail << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";

        std::vector< double > vec_gmpm;
        std::vector< double > vec_gm  ;
        std::vector< double > vec_pm  ;
        std::vector< double > vec_tail;

        std::vector< std::pair< double, std::string >> temp_gmpm;
        std::vector< std::pair< double, std::string >> temp_gm  ;
        std::vector< std::pair< double, std::string >> temp_pm  ;
        std::vector< std::pair< double, std::string >> temp_tail;

        std::tuple< double, std::size_t, std::size_t > df_tuple_gmpm;
        std::tuple< double, std::size_t, std::size_t > df_tuple_gm  ;
        std::tuple< double, std::size_t, std::size_t > df_tuple_pm  ;
        std::tuple< double, std::size_t, std::size_t > df_tuple_tail;

        double gm = 0.0;
        double pm = 0.0;
        double sum_loading = 0.0;

        std::string res_gmpm;
        std::string res_gm  ;
        std::string res_pm  ;
        std::string res_tail;

        for( auto& anno : ano_len_idx.first )
        {
            vec_gmpm.clear();
            vec_gm  .clear();
            vec_pm  .clear();
            vec_tail.clear();

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
                vec_gm  .emplace_back( gm );
                vec_pm  .emplace_back( pm );
                vec_tail.emplace_back( pm < 1 ? 0 : ( pm * 100 / ( gm + pm )));
            }

            sum_loading = get_sum( vec_gmpm );

            df_tuple_gmpm = get_difference( vec_gmpm );
            df_tuple_gm   = get_difference( vec_gm   );
            df_tuple_pm   = get_difference( vec_pm   );
            df_tuple_tail = get_difference( vec_tail );

            res_gmpm = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_gm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_pm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tail = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_gmpm += std::to_string( std::get<0>( df_tuple_gmpm )) + "\t" + bed_samples[ std::get<1>( df_tuple_gmpm )].first + ":" + bed_samples[ std::get<2>( df_tuple_gmpm )].first;
            res_gm   += std::to_string( std::get<0>( df_tuple_gm   )) + "\t" + bed_samples[ std::get<1>( df_tuple_gm   )].first + ":" + bed_samples[ std::get<2>( df_tuple_gm   )].first;
            res_pm   += std::to_string( std::get<0>( df_tuple_pm   )) + "\t" + bed_samples[ std::get<1>( df_tuple_pm   )].first + ":" + bed_samples[ std::get<2>( df_tuple_pm   )].first;
            res_tail += std::to_string( std::get<0>( df_tuple_tail )) + "\t" + bed_samples[ std::get<1>( df_tuple_tail )].first + ":" + bed_samples[ std::get<2>( df_tuple_tail )].first;

            temp_gmpm.emplace_back( std::make_pair( sum_loading, res_gmpm ));
            temp_gm  .emplace_back( std::make_pair( sum_loading, res_gm   ));
            temp_pm  .emplace_back( std::make_pair( sum_loading, res_pm   ));
            temp_tail.emplace_back( std::make_pair( sum_loading, res_tail ));
        }

        sort_difference( temp_gmpm );
        sort_difference( temp_gm   );
        sort_difference( temp_pm   );
        sort_difference( temp_tail );

        for( auto& output : temp_gmpm ) out_gmpm << output.second;
        for( auto& output : temp_gm   ) out_gm   << output.second;
        for( auto& output : temp_pm   ) out_pm   << output.second;
        for( auto& output : temp_tail ) out_tail << output.second;

        out_gmpm << "\n";
        out_gm   << "\n";
        out_pm   << "\n";
        out_tail << "\n";

        out_gmpm.close();
        out_gm  .close();
        out_pm  .close();
        out_tail.close();
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

            for( auto& lens : len_vec[ smp ] ) if( max < lens.second )
            {
                len = lens.first;
                max = lens.second;
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

            if( max.first < len_vec[ smp ] ) max = std::make_pair( len_vec[ smp ], smp );
            if( min.first > len_vec[ smp ] ) min = std::make_pair( len_vec[ smp ], smp );
        }

        return { max.first - min.first, max.second, min.second };
    }

    void output_length_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        std::ofstream out_gmpm( output_path + "LengthDifference_GMPM.tsv"    );
        std::ofstream out_gm  ( output_path + "LengthDifference_GM.tsv"      );
        std::ofstream out_pm  ( output_path + "LengthDifference_PM.tsv"      );
        std::ofstream out_tail( output_path + "LengthDifference_Tailing.tsv" );

        out_gmpm << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";
        out_gm   << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";
        out_pm   << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";
        out_tail << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::vector< std::size_t > len_gmpm;
        std::vector< std::size_t > len_gm  ;
        std::vector< std::size_t > len_pm  ;
        std::vector< std::size_t > len_tail;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::size_t, std::tuple< double, double >>>> smp_anno_len_gmpms;

        std::vector< std::map< std::size_t, double >> vec_len_gmpm( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_gm  ( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_pm  ( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_tail( bed_samples.size(), std::map< std::size_t, double >() );

        std::vector< std::pair< double, std::string >> temp_gmpm;
        std::vector< std::pair< double, std::string >> temp_gm  ;
        std::vector< std::pair< double, std::string >> temp_pm  ;
        std::vector< std::pair< double, std::string >> temp_tail;

        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_gmpm;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_gm  ;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_pm  ;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_tail;

        std::string res_gmpm;
        std::string res_gm  ;
        std::string res_pm  ;
        std::string res_tail;

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
            vec_len_gmpm = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );
            vec_len_gm   = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );
            vec_len_pm   = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );
            vec_len_tail = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& len : smp_anno_len_gmpms[ smp ][ anno ] )
                {
                    if( len.first == 0 )
                    {
                        vec_loading.emplace_back( std::get<0>( len.second ) );
                        continue;
                    }

                    vec_len_gmpm[ smp ][ len.first ] = std::get<0>( len.second ) + std::get<1>( len.second );
                    vec_len_gm  [ smp ][ len.first ] = std::get<0>( len.second );
                    vec_len_pm  [ smp ][ len.first ] = std::get<1>( len.second );
                    vec_len_tail[ smp ][ len.first ] =
                        ( std::get<1>( len.second ) < 1 ? 0 :( std::get<1>( len.second ) * 100 /( std::get<0>( len.second ) + std::get<1>( len.second ))));
                }
            }

            sum_loading = get_sum( vec_loading );

            len_gmpm = get_lengths( vec_len_gmpm );
            len_gm   = get_lengths( vec_len_gm   );
            len_pm   = get_lengths( vec_len_pm   );
            len_tail = get_lengths( vec_len_tail );

            df_tuple_gmpm = get_length_difference( len_gmpm );
            df_tuple_gm   = get_length_difference( len_gm   );
            df_tuple_pm   = get_length_difference( len_pm   );
            df_tuple_tail = get_length_difference( len_tail );

            res_gmpm = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_gm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_pm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tail = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_gmpm += std::to_string( std::get<0>( df_tuple_gmpm )) + "\t";
            res_gm   += std::to_string( std::get<0>( df_tuple_gm   )) + "\t";
            res_pm   += std::to_string( std::get<0>( df_tuple_pm   )) + "\t";
            res_tail += std::to_string( std::get<0>( df_tuple_tail )) + "\t";

            res_gmpm += bed_samples[ std::get<1>( df_tuple_gmpm )].first + ":" + bed_samples[ std::get<2>( df_tuple_gmpm )].first + "\t";
            res_gm   += bed_samples[ std::get<1>( df_tuple_gm   )].first + ":" + bed_samples[ std::get<2>( df_tuple_gm   )].first + "\t";
            res_pm   += bed_samples[ std::get<1>( df_tuple_pm   )].first + ":" + bed_samples[ std::get<2>( df_tuple_pm   )].first + "\t";
            res_tail += bed_samples[ std::get<1>( df_tuple_tail )].first + ":" + bed_samples[ std::get<2>( df_tuple_tail )].first + "\t";

            res_gmpm += std::to_string( len_gmpm[ std::get<1>( df_tuple_gmpm )]) + ":" + std::to_string( len_gmpm[ std::get<2>( df_tuple_gmpm )]);
            res_gm   += std::to_string( len_gm  [ std::get<1>( df_tuple_gm   )]) + ":" + std::to_string( len_gm  [ std::get<2>( df_tuple_gm   )]);
            res_pm   += std::to_string( len_pm  [ std::get<1>( df_tuple_pm   )]) + ":" + std::to_string( len_pm  [ std::get<2>( df_tuple_pm   )]);
            res_tail += std::to_string( len_tail[ std::get<1>( df_tuple_tail )]) + ":" + std::to_string( len_tail[ std::get<2>( df_tuple_tail )]);

            temp_gmpm.emplace_back( std::make_pair( sum_loading, res_gmpm ));
            temp_gm  .emplace_back( std::make_pair( sum_loading, res_gm   ));
            temp_pm  .emplace_back( std::make_pair( sum_loading, res_pm   ));
            temp_tail.emplace_back( std::make_pair( sum_loading, res_tail ));
        }

        sort_difference( temp_gmpm );
        sort_difference( temp_gm   );
        sort_difference( temp_pm   );
        sort_difference( temp_tail );

        for( auto& output : temp_gmpm ) out_gmpm << output.second;
        for( auto& output : temp_gm   ) out_gm   << output.second;
        for( auto& output : temp_pm   ) out_pm   << output.second;
        for( auto& output : temp_tail ) out_tail << output.second;

        out_gmpm << "\n";
        out_gm   << "\n";
        out_pm   << "\n";
        out_tail << "\n";

        out_gmpm.close();
        out_gm  .close();
        out_pm  .close();
        out_tail.close();
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

    void output_arms_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        std::ofstream out_gmpm( output_path + "ArmDifference_GMPM.tsv"    );
        std::ofstream out_gm  ( output_path + "ArmDifference_GM.tsv"      );
        std::ofstream out_pm  ( output_path + "ArmDifference_PM.tsv"      );
        std::ofstream out_tail( output_path + "ArmDifference_Tailing.tsv" );

        out_gmpm << "Annotation\tTotal\tisLoadingDifference";
        out_gm   << "Annotation\tTotal\tisLoadingDifference";
        out_pm   << "Annotation\tTotal\tisLoadingDifference";
        out_tail << "Annotation\tTotal\tisLoadingDifference";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_gmpm << "\tSample1";
            else out_gmpm << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_gm   << "\tSample1";
            else out_gm   << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_pm   << "\tSample1";
            else out_pm   << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_tail << "\tArm1";
            else out_tail << ":Arm" << smp + 1 ;

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::string, std::tuple< double, double >>>> smp_anno_arm_gmpms;

        std::vector< std::map< std::string, double >> vec_arm_gmpm( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::map< std::string, double >> vec_arm_gm  ( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::map< std::string, double >> vec_arm_pm  ( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::map< std::string, double >> vec_arm_tail( bed_samples.size(), std::map< std::string, double >() );

        std::vector< std::pair< std::string, double >> arm_gmpm;
        std::vector< std::pair< std::string, double >> arm_gm  ;
        std::vector< std::pair< std::string, double >> arm_pm  ;
        std::vector< std::pair< std::string, double >> arm_tail;

        std::vector< std::pair< double, std::string >> temp_gmpm;
        std::vector< std::pair< double, std::string >> temp_gm  ;
        std::vector< std::pair< double, std::string >> temp_pm  ;
        std::vector< std::pair< double, std::string >> temp_tail;

        std::string is_diff_gmpm;
        std::string is_diff_gm  ;
        std::string is_diff_pm  ;
        std::string is_diff_tail;

        std::string res_gmpm;
        std::string res_gm  ;
        std::string res_pm  ;
        std::string res_tail;

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

            vec_arm_gmpm = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );
            vec_arm_gm   = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );
            vec_arm_pm   = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );
            vec_arm_tail = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& arm : smp_anno_arm_gmpms[ smp ][ anno ] )
                {
                    vec_loading.emplace_back( std::get<0>( arm.second ) + std::get<1>( arm.second ));

                    vec_arm_gmpm[ smp ][ arm.first ] =   std::get<0>( arm.second ) + std::get<1>( arm.second );
                    vec_arm_gm  [ smp ][ arm.first ] =   std::get<0>( arm.second ); 
                    vec_arm_pm  [ smp ][ arm.first ] =   std::get<1>( arm.second );
                    vec_arm_tail[ smp ][ arm.first ] = ( std::get<1>( arm.second ) < 1 ? 0 :( std::get<1>( arm.second ) * 100 /( std::get<0>( arm.second ) + std::get<1>( arm.second ))));
                }
            }

            sum_loading = get_sum(  vec_loading );

            arm_gmpm = get_arms( vec_arm_gmpm );
            arm_gm   = get_arms( vec_arm_gm   );
            arm_pm   = get_arms( vec_arm_pm   );
            arm_tail = get_arms( vec_arm_tail );

            is_diff_gmpm = get_arm_difference( arm_gmpm );
            is_diff_gm   = get_arm_difference( arm_gm   );
            is_diff_pm   = get_arm_difference( arm_pm   );
            is_diff_tail = get_arm_difference( arm_tail );

            res_gmpm = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_gmpm;
            res_gm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_gm  ;
            res_pm   = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_pm  ;
            res_tail = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_tail;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_gmpm += "\t" + bed_samples[ smp ].first;
                else res_gmpm += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_gmpm += "\t" + arm_gmpm[ smp ].first;
                else res_gmpm += ":" + arm_gmpm[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_gm   += "\t" + bed_samples[ smp ].first;
                else res_gm   += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_gm   += "\t" + arm_gm  [ smp ].first;
                else res_gm   += ":" + arm_gm  [ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_pm   += "\t" + bed_samples[ smp ].first;
                else res_pm   += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_pm   += "\t" + arm_pm  [ smp ].first;
                else res_pm   += ":" + arm_pm  [ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tail += "\t" + bed_samples[ smp ].first;
                else res_tail += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tail += "\t" + arm_tail[ smp ].first;
                else res_tail += ":" + arm_tail[ smp ].first;

            temp_gmpm.emplace_back( std::make_pair( sum_loading, res_gmpm ));
            temp_gm  .emplace_back( std::make_pair( sum_loading, res_gm   ));
            temp_pm  .emplace_back( std::make_pair( sum_loading, res_pm   ));
            temp_tail.emplace_back( std::make_pair( sum_loading, res_tail ));
        }

        sort_difference( temp_gmpm );
        sort_difference( temp_gm   );
        sort_difference( temp_pm   );
        sort_difference( temp_tail );

        for( auto& output : temp_gmpm ) out_gmpm << output.second;
        for( auto& output : temp_gm   ) out_gm   << output.second;
        for( auto& output : temp_pm   ) out_pm   << output.second;
        for( auto& output : temp_tail ) out_tail << output.second;

        out_gmpm << "\n";
        out_gm   << "\n";
        out_pm   << "\n";
        out_tail << "\n";

        out_gmpm.close();
        out_gm  .close();
        out_pm  .close();
        out_tail.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
