#pragma once

namespace ago {
namespace algorithm {

class MetaAnalyzer
{
  public:

    MetaAnalyzer()
    {}

    void quantile_transfer(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
        , std::vector< QuantileDataType<> >& qres
        , const std::vector< std::string >& index
        , const char* type
        , const char* gmpm
    )
    {
        std::string read_count = std::string( type ) + "_" + std::string( gmpm ) + "_read_count";
        std::string ppm        = std::string( type ) + "_" + std::string( gmpm ) + "_ppm";

        ago::engine::DataPool::anno_len_samples rd_res;
        ago::engine::DataPool::anno_len_samples pp_res;

        size_t sample_count = 0;

        // if( type == ".LenDist" )
        // {
        //     std::ofstream output( "test.log" );
        //     for( auto& samp : analyzer_result_samples )
        //         for( auto& anno_len : samp.second )
        //             for( auto& anno : anno_len )
        //             {
        //                 for( auto& len : anno.second )
        //                     output << samp.first << "\t" << anno.first << "\t" << len.first << "\t" << len.second << "\n";
        //                 output << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        //             }
        //     std::ofstream output( "test.log" );
        //     output.close();
        //     output.close();
        // }
        
        std::map< std::size_t, std::vector< int >> lens_index;
        
        for( std::size_t i = 0; i < analyzer_result_samples.size(); ++i )
        {
            for( auto& annos : analyzer_result_samples[i].second )
            {
                if( annos.find( read_count ) != annos.end() )
                {
                    std::size_t idx = 0;
                    for( auto& len : annos.begin()->second )
                    {
                        if( lens_index.find( len.second ) == lens_index.end() )
                            lens_index[ len.second ] = std::vector< int >( analyzer_result_samples.size(), -1 );

                        lens_index[ len.second ][i] = idx; 
                        idx++;
                    }
                }
            }
        }

        // for( auto& len : lens_index )
        // {
        //     std::cerr << len.first;
        //     for( auto& smp : len.second )
        //         std::cerr << "\t" << smp;
        //     std::cerr << "\n";
        // }

        std::vector< double > head_vec;
        for( auto& len : lens_index )
        {
            if( len.first == 0 ) continue;
            head_vec.emplace_back( len.first );
        }
        head_vec.emplace_back( 0 );

        for( std::size_t smp = 0; smp < analyzer_result_samples.size(); ++smp )
        {
            auto& sample = analyzer_result_samples[ smp ];

            std::map< std::string, std::vector< double >> rd_anno_map;
            std::map< std::string, std::vector< double >> pp_anno_map;

            for( auto& annos : sample.second )
            {
                if( annos.find( read_count ) != annos.end() )
                {
                    rd_anno_map.emplace( read_count, head_vec );
                    pp_anno_map.emplace( ppm, head_vec );

                    size_t anno_count = 0;

                    for( auto& idx : index )
                    {
                        auto it = annos.find( idx );

                        if( it != annos.end() )
                        {
                            double average = 0;
                            double ppm_sum = 0;

                            if( gmpm == "GMPM" )
                            {
                                average = qres[ sample_count ].value_[ anno_count ].first / it->second.find( "SUM_LEN" )->second;
                            }
                            else
                            {
                                for( auto& ann : sample.second )
                                {
                                    if( ann.find( std::string( type ) + "_GMPM_read_count" ) != ann.end() )
                                    {
                                         ppm_sum = ann.find( idx )->second.find( "SUM_LEN" )->second;
                                    }
                                }
                                average = qres[ sample_count ].value_[ anno_count ].first / ppm_sum;
                            }

                            std::vector< double > rd_vec = std::vector< double >( head_vec.size(), 0 );
                            std::vector< double > pp_vec = std::vector< double >( head_vec.size(), 0 );

                            int lenidx = 0;
                            int shift = 0;
                            for( auto& len : it->second )
                            {
                                if( lens_index[ head_vec[ lenidx ]][ smp ] == -1 ) shift++;
                                rd_vec[ lens_index[ head_vec[ lenidx + shift ]][ smp ] + shift ] = len.second;
                                pp_vec[ lens_index[ head_vec[ lenidx + shift ]][ smp ] + shift ] = len.second * average;

                                // if( type == ".LenDist" && idx == "miRNA" )
                                //     std::cerr << len.first << "\t"
                                //               << len.second << "\t"
                                //               << lenidx << "+" << shift << "\t"
                                //               << head_vec[ lenidx + shift ] << "\t" 
                                //               << lens_index[ head_vec[ lenidx + shift ]][ smp ] << "+" << shift << "\n";
                                lenidx++;
                            }

                            rd_anno_map.emplace( idx, rd_vec );
                            pp_anno_map.emplace( idx, pp_vec );
                        }
                        else
                        {
                            rd_anno_map.emplace( idx, std::vector< double >( head_vec.size(), 0 ));
                            pp_anno_map.emplace( idx, std::vector< double >( head_vec.size(), 0 ));
                        }

                        anno_count++;
                    }
                }
            }

            std::string anno_sum;

            if( type == ".LenDist" )
            {
                anno_sum = "SUM_ANNO";
            }
            else
            {
                anno_sum = "SUM_ANNO_SUM_ANNO";
            }

            rd_anno_map.emplace( anno_sum, get_sum_vec( rd_anno_map ));
            pp_anno_map.emplace( anno_sum, get_sum_vec( pp_anno_map ));

            // emplace_lens( rd_anno_map );
            // emplace_lens( pp_anno_map );

            rd_res.emplace( sample.first, rd_anno_map );
            pp_res.emplace( sample.first, pp_anno_map );

            sample_count++;
        }

        all_results.emplace_back( read_count, rd_res );
        all_results.emplace_back( ppm, pp_res );

        // if( type == ".LenDist" )
        // {
        //     std::ofstream output( "test.log" );
        //     for( auto& result_type : all_results )
        //         for( auto& samp : result_type.second )
        //             for( auto& anno_len : samp.second )
        //                 for( auto& len : anno_len.second )
        //                     output << result_type.first << "\t" << samp.first << "\t" << anno_len.first << "\t" << len << "\n";
        //     output.close();
        // }
    }

    std::vector< double > get_sum_vec( const std::map< std::string, std::vector< double >>& annos )
    {
        std::vector< double > sum_vec;

        for( int i = 0; i < annos.begin()->second.size(); ++i )
        {
            double sum = 0;
            int flag = 0;

            for( auto& anno : annos )
            {
                if( flag != 0 && !anno.second.empty() )
                {
                    sum += anno.second[i];
                }

                flag++;
            }

            sum_vec.push_back( sum );
        }

        return sum_vec;
    }

    void emplace_lens( std::map< std::string, std::vector< double >>& annos )
    {
        std::vector< double > len_index{
            15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 0
        //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 16
        };

        std::map< size_t, size_t > len_emplace;

        int j = -1;
        size_t k = 0;

        for( size_t i = 0; i < len_index.size(); ++i )
        {
            j++;

            if( len_index[ i ] != annos.begin()->second[ j ] )
            {
                if( annos.begin()->second[ j ] != 0 )
                {
                    len_emplace.emplace( j, len_index[ i ] );
                }
                else
                {
                    switch( k )
                    {
                        case 0: k = j; break;
                        default: k++;
                    }

                    len_emplace.emplace( k, len_index[ i ] );
                }

                j--;
            }
        }

        size_t flag = 0;
        for( auto& len : annos )
        {
            if( len.second.empty() )
                continue;

            size_t shift = 0;
            size_t size = len.second.size() - 1;
            size_t size_shift = 0;

            for( auto& emplace : len_emplace )
            {
                if( emplace.first >= size )
                {
                    switch( flag )
                    {
                        case 0:
                            len.second.emplace( len.second.begin() + emplace.first + size_shift, emplace.second );
                            break;

                        default:
                            len.second.emplace( len.second.begin() + emplace.first + size_shift, 0 );
                    }
                }
                else
                {
                    size_shift++;
                    switch( flag )
                    {
                        case 0:
                            len.second.emplace( len.second.begin() + emplace.first + shift, emplace.second );
                            break;

                        default:
                            len.second.emplace( len.second.begin() + emplace.first + shift, 0 );
                    }

                    if( shift < 2 ) shift++;
                }
            }

            flag++;
        }
    } 

    void tailing_ratio(
            std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
          , std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
          , const char* type
          )
    {
        ago::engine::DataPool::anno_len_samples res;

        std::map< std::string, size_t > anno_index;
        std::map< size_t, size_t > len_index;

        for( auto& sample : analyzer_result_samples )
        {
            for( auto& annos : sample.second )
            {
                if( annos.find( type ) != annos.end() )
                {
                    int flag = 0;

                    for( auto& pair : annos )
                    {
                        switch( flag )
                        {
                            case 0:
                                for( auto& length : pair.second )
                                {
                                    len_index.emplace( length.second, 0 );
                                }
                                break;

                            default:
                                anno_index.emplace( pair.first, 0 );
                        }

                        flag++;
                    }
                }
            }
        }

        for( auto& sample : analyzer_result_samples )
        {
            for( auto& annos : sample.second )
            {
                if( annos.find( type ) != annos.end() )
                {
                    std::map< std::string, std::vector< double >> ratio_table;
                    {
                        std::vector< double > ratios;

                        for( auto& len : len_index )
                        {
                            if( len.first != 0 )
                            {
                                ratios.push_back( len.first );
                            }
                        }

                        ratios.push_back( len_index.begin()->first );
                        ratio_table.emplace( type, ratios );
                    }

                    for( auto& anno : anno_index )
                    {
                        std::vector< double > ratios;

                        std::map< std::string, double >::iterator len_it;
                        std::map< std::string, std::map< std::string, double >>::iterator ratio_it;

                        ratio_it = annos.find( anno.first );

                        if( ratio_it != annos.end() )
                        {
                            for( auto& len : len_index )
                            {
                                if( len.first != 0 )
                                {
                                    len_it = ratio_it->second.find( std::to_string( len.first ));

                                    if( len_it != ratio_it->second.end() )
                                    {
                                        ratios.push_back( len_it->second );
                                    }
                                    else
                                    {
                                        ratios.push_back( len.second );
                                    }
                                }
                            }

                            len_it = ratio_it->second.find( "SUM_LEN" );
                            ratios.push_back( len_it->second );
                        }
                        else
                        {
                            for( auto& len : len_index )
                                ratios.push_back(0);
                        }

                        ratio_table.emplace( anno.first, ratios );
                    }

                    res.emplace( sample.first, ratio_table );
                }
            }
        }
        all_results.emplace_back( type, res );
    }

    void tailing_ratio_for_pm(
            std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
          , std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
          , const char* type
          , const ago::engine::DataPool::anno_len_samples& lendist_gm_readcount
          , const ago::engine::DataPool::anno_len_samples& lendist_gm_ppm
          , const ago::engine::DataPool::anno_len_samples& lendist_pm_readcount
          , const ago::engine::DataPool::anno_len_samples& lendist_pm_ppm
    )
    {
        ago::engine::DataPool::anno_len_samples res;

		std::map< std::string, std::map< std::string, std::vector< double >>> gm_pm_tail_readcount;
		std::map< std::string, std::map< std::string, std::vector< double >>> gm_pm_tail_ppm ;

        set_lendist_pm( analyzer_result_samples, type, gm_pm_tail_readcount, lendist_pm_readcount, ".LenDist_PM_read_count", ".GMPMTR_read_count" );
        set_lendist_gm( gm_pm_tail_readcount, lendist_gm_readcount, ".LenDist_GM_read_count" );
        set_lendist_pm( analyzer_result_samples, type, gm_pm_tail_ppm, lendist_pm_ppm, ".LenDist_PM_ppm", ".GMPMTR_ppm" );
        set_lendist_gm( gm_pm_tail_ppm, lendist_gm_ppm, ".LenDist_GM_ppm" );

        all_results.emplace_back( ".GMPMTR_read_count", gm_pm_tail_readcount );
        all_results.emplace_back( ".GMPMTR_ppm", gm_pm_tail_ppm );
    }

    void set_lendist_pm(
          std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
        , const char* type
        , std::map< std::string, std::map< std::string, std::vector< double >>>& gmpm_tail_rdpp
        , const ago::engine::DataPool::anno_len_samples& lendist_pm
        , const char* lendist_pm_type
        , const char* gmpmtr_type
    )
    {
        for( auto& pm_value : lendist_pm )
        {
            std::map< std::string, std::vector< double >> tail_length;
            std::vector< double > pmlen_index;

            for( auto& anno : pm_value.second )
            {
                if( anno.first == lendist_pm_type )
                {
                    //                          A  C  G  T  X(other)
                    std::vector< double > atcg{ 0, 1, 2, 3, 4 };
                    tail_length.emplace( gmpmtr_type, atcg );

                    for( auto& len : anno.second )
                        pmlen_index.push_back( len );
                }

                if( anno.first == "miRNA" )
                {
                    for( auto& sample : analyzer_result_samples )
                    {
                        for( auto& annos : sample.second )
                        {
                            if( annos.find( type ) != annos.end() )
                            {
                                for( int len = 0; len < anno.second.size(); len++ )
                                {
                                    if( (int)pmlen_index[ len ] == 0 )
                                        continue;

                                    auto anno_it = annos.find( std::to_string( (int)pmlen_index[ len ] ));

                                    if( anno_it != annos.end() )
                                    {
                                        std::vector< double > atcg;

                                        for( auto& n_tail : anno_it->second )
                                            atcg.push_back( anno.second[ len ] * n_tail.second );

                                        if( anno_it->second.size() == 4 )
                                            atcg.push_back( 0 );

                                        tail_length.emplace( std::to_string( (int)pmlen_index[ len ] ), atcg );
                                    }
                                    else
                                    {
                                        tail_length.emplace( std::to_string( (int)pmlen_index[ len ] ), std::vector<double>( 5, 0 ));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            gmpm_tail_rdpp.emplace( pm_value.first, tail_length);
        }
    }

    void set_lendist_gm(
          std::map< std::string, std::map< std::string, std::vector< double >>>& gmpm_tail_rdpp
        , const ago::engine::DataPool::anno_len_samples& lendist_gm
        , const char* lendist_gm_type
    )
    {
        for( auto& gm_value : lendist_gm )
        {
            std::vector< double > gmlen_index;

            for( auto& anno : gm_value.second )
            {
                if( anno.first == lendist_gm_type )
                {
                    for( auto& len : anno.second )
                    {
                        if( len != 0 )
                        {
                            gmlen_index.push_back( len );
                        }
                    }
                    //                                                                      GM
                    gmpm_tail_rdpp.find( gm_value.first )->second.begin()->second.push_back( 5 );
                }

                if( anno.first == "miRNA" )
                {
                    for( int len = 0; len < anno.second.size(); len++ )
                    {
                        if( (int)gmlen_index[ len ] == 0)
                            continue;

                        auto gmpm_it = gmpm_tail_rdpp.find( gm_value.first )->second.find( std::to_string( (int)gmlen_index[ len ] ));

                        if( gmpm_it != gmpm_tail_rdpp.find( gm_value.first )->second.end() )
                        {
                            gmpm_it->second.push_back(anno.second[len]);
                        }
                    }
                }
            }
        }
    }

    void tailing_ratio_for_mir(
            std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
          , std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
          , const char* type
          , const ago::engine::DataPool::anno_len_samples& mirdist_gmpm_ppm
    )
    {
        ago::engine::DataPool::anno_len_samples res;

        for( auto& sample : analyzer_result_samples )
        {
            for( auto& annos : sample.second )
            {
                if( annos.find( type ) != annos.end() )
                {
                    std::map< std::string, std::vector< double >> anno_value;

                    //                                  other
                    //                          A  C  G  T  X  Z
                    //                                         GM
                    std::vector< double > ATGC{ 0, 1, 2, 3, 4, 5};
                    anno_value.emplace( ".MirTail_ppm", ATGC );

                    for( auto& anno : annos )
                    {
                        std::vector< std::string > anno_len;
                        boost::iter_split( anno_len, anno.first, boost::algorithm::first_finder( ":" ));

                        for( auto& mir_sample : mirdist_gmpm_ppm )
                        {
                            if( mir_sample.first == sample.first &&
                                    mir_sample.second.find( anno_len[0] ) != mir_sample.second.end() )
                            {
                                for( int l = 0; l < mir_sample.second.find( ".MirDist_GMPM_ppm" )->second.size(); ++l )
                                {
                                    if( mir_sample.second.find( ".MirDist_GMPM_ppm" )->second[l] == std::stod( anno_len[1] ))
                                    {
                                        std::vector< double > value;

                                        for( auto& val : anno.second )
                                        {
                                            if( mir_sample.second.find( anno_len[0] )->second.empty() )
                                                continue;

                                            val.second = val.second * mir_sample.second.find( anno_len[0] )-> second[l];
                                            value.push_back( val.second );
                                        }

                                        anno_value.emplace( anno.first, value );
                                    }
                                }
                            }
                        }
                    }
                    res.emplace( sample.first, anno_value );
                }
            }
        }
        all_results.emplace_back( ".MirTail_ppm", res );
    }

};

} // end of namespace algorithm
} // end of namespace ago
