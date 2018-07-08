#pragma once
#include <cmath>
#include <boost/math/special_functions/beta.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDifferential
{
  public:

    GeneTypeAnalyzerDifferential()
    {}

    static std::vector< std::vector< std::pair< std::size_t, std::string >>>
        get_smp_compare( const std::vector< BedSampleType >& bed_samples )
    {
        std::vector< std::pair< std::size_t, std::string >> compares_temp;
        std::vector< std::vector< std::pair< std::size_t, std::string >>> compares;

        std::set< std::pair< std::size_t, std::string >> smp_check_temp;
        std::set< std::set< std::pair< std::size_t, std::string >>> smp_check;

        for( std::size_t s1 = 0; s1 < bed_samples.size(); ++s1 )
        {
            for( std::size_t s2 = 0; s2 < bed_samples.size(); ++s2 )
            {
                if( s1 == s2 ) continue;
                smp_check_temp.emplace( std::make_pair( s1, bed_samples[s1].first ));
                smp_check_temp.emplace( std::make_pair( s2, bed_samples[s2].first ));
                smp_check.emplace( smp_check_temp );
                smp_check_temp.clear();
            }
        }

        for( auto& compare : smp_check )
        {
            for( auto& smp : compare ) compares_temp.emplace_back( smp );
            compares.emplace_back( compares_temp );
            compares_temp.clear();
        }

        return compares;
    }

    static std::string output_header( const std::string& s1, const std::string& s2 )
    {
        return "Annotation\tTotal\t"
            + s1 + ":" + s2 + ":GMPM\t"
            + s1 + ":" + s2 + ":GM\t"
            + s1 + ":" + s2 + ":PM\t"
            + s1 + ":" + s2 + ":Atail\t"
            + s1 + ":" + s2 + ":Ctail\t"
            + s1 + ":" + s2 + ":Gtail\t"
            + s1 + ":" + s2 + ":Ttail\t"
            + s1 + ":" + s2 + ":Other\t"
            + s2 + ":" + s1 + ":GMPM\t"
            + s2 + ":" + s1 + ":GM\t"
            + s2 + ":" + s1 + ":PM\t"
            + s2 + ":" + s1 + ":Atail\t"
            + s2 + ":" + s1 + ":Ctail\t"
            + s2 + ":" + s1 + ":Gtail\t"
            + s2 + ":" + s1 + ":Ttail\t"
            + s2 + ":" + s1 + ":Other";
    }

    static double get_value( auto& anno_table_tail, auto& idx, auto& anno, auto& lens, const auto& token )
    {
        double res = 0.0;

        for( std::size_t i = 0; i < 5; i++ )
        {
            switch( token )
            {
                case 'A' : i = 0; break;
                case 'C' : i = 1; break;
                case 'T' : i = 2; break;
                case 'G' : i = 3; break;
                case 'O' : i = 4; break;
                case 'M' : i = 5; break;
            }

            if( anno_table_tail[ idx ][i].find( anno ) != anno_table_tail[ idx ][i].end() )
                for( auto& len : lens )
                    if( anno_table_tail[ idx ][i][ anno ].find( len ) != anno_table_tail[ idx ][i][ anno ].end() )
                        res += anno_table_tail[ idx ][i][ anno ][ len ];

            if( token != 'P' ) break;
        }

        return res;
    }

    /*
     * P value = 1 - CDF
     * CDF = CDF( X, A, B )
     * X = Fold Change / ( 1 + Fold Change )
     * A = 3.26 * ( Log( Count1 ) / CountMax )
     * B = 1.63 * ( Log( Count2 ) / CountMax )
     * CountMax = Max( Log( Count1 ), Log( Count2 ))
     */
    static std::pair< double, double > get_differential( const auto& s1, const auto& s2 )
    {
        double fold_change = s1 / s2;
        double count1 = std::log( s1 );
        double count2 = std::log( s2 );
        double count_max = ( count1 >= count2 ? count1 : count2 );
        double p_value = 1 - boost::math::ibeta(( 3.26 * ( count1 / count_max )), ( 1.63 * ( count2 / count_max )), ( fold_change / ( 1 + fold_change )));
        return std::make_pair( fold_change, p_value );
    }

    static double get_smallest_pvalue( auto& diff_vec )
    {
        double p_value = 999;

        for( auto& diff : diff_vec )
            if( diff.second < p_value ) p_value = diff.second;

        return p_value;
    }

    static std::string output_formation( auto& diff_vec )
    {
        std::string out = "";

        for( auto& diff : diff_vec )
            out += "\t" + std::to_string( diff.first ) + ":" + std::to_string( diff.second );

        return out;
    }

    static void sort_differential( std::vector< std::tuple< double, double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::tuple< double, double, std::string >& a, const std::tuple< double, double, std::string >& b )
            {
                     if( std::get<0>(a) != std::get<0>(b) ) return std::get<0>(a) > std::get<0>(b);
                else if( std::get<1>(a) != std::get<1>(b) ) return std::get<1>(a) < std::get<1>(b);
                else                                        return std::get<2>(a) > std::get<2>(b);
            });
    }

    static void output_loading_differential(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        std::ofstream output;
        for( auto& compare : get_smp_compare( bed_samples ))
        {
            auto& s1 = compare[0];
            auto& s2 = compare[1];

            output.open( output_path + "LoadingDifferential_" + s1.second + "_" + s2.second + ".text" );
            output << output_header( s1.second, s2.second );

            double gm1, gm2, pm1, pm2, at1, at2, ct1, ct2, gt1, gt2, tt1, tt2, ot1, ot2;
            std::vector< std::pair< double, double >> diff_vec;
            std::vector< std::tuple< double, double, std::string >> out_temp;
            //                       total  smallestPvalue

            for( auto& anno : ano_len_idx.first )
            {
                gm1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'M' );
                gm2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'M' );
                pm1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'P' );
                pm2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'P' );
                at1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'A' );
                at2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'A' );
                ct1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'C' );
                ct2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'C' );
                gt1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'G' );
                gt2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'G' );
                tt1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'T' );
                tt2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'T' );
                ot1 = get_value( anno_table_tail, s1.first, anno, ano_len_idx.second, 'O' );
                ot2 = get_value( anno_table_tail, s2.first, anno, ano_len_idx.second, 'O' );

                diff_vec.emplace_back( get_differential( gm1 + pm1 + 1, gm2 + pm2  + 1 ));
                diff_vec.emplace_back( get_differential( gm1       + 1, gm2        + 1 ));
                diff_vec.emplace_back( get_differential(       pm1 + 1,       pm2  + 1 ));
                diff_vec.emplace_back( get_differential(    at1    + 1,    at2     + 1 ));
                diff_vec.emplace_back( get_differential(    ct1    + 1,    ct2     + 1 ));
                diff_vec.emplace_back( get_differential(    gt1    + 1,    gt2     + 1 ));
                diff_vec.emplace_back( get_differential(    tt1    + 1,    tt2     + 1 ));
                diff_vec.emplace_back( get_differential(    ot1    + 1,    ot2     + 1 ));

                diff_vec.emplace_back( get_differential( gm2 + pm2 + 1, gm1 + pm1  + 1 ));
                diff_vec.emplace_back( get_differential( gm2       + 1, gm1        + 1 ));
                diff_vec.emplace_back( get_differential(       pm2 + 1,       pm1  + 1 ));
                diff_vec.emplace_back( get_differential(    at2    + 1,    at1     + 1 ));
                diff_vec.emplace_back( get_differential(    ct2    + 1,    ct1     + 1 ));
                diff_vec.emplace_back( get_differential(    gt2    + 1,    gt1     + 1 ));
                diff_vec.emplace_back( get_differential(    tt2    + 1,    tt1     + 1 ));
                diff_vec.emplace_back( get_differential(    ot2    + 1,    ot1     + 1 ));

                out_temp.emplace_back( std::make_tuple(
                        gm1 + pm1 + gm2 + pm2,
                        get_smallest_pvalue( diff_vec ),
                        "\n" + anno + "\t"
                             + std::to_string( gm1 + pm1 + gm2 + pm2 )
                             + output_formation( diff_vec )
                        ));

                diff_vec.clear();
            }

            sort_differential( out_temp );
            for( auto& temp : out_temp )
                output << std::get<2>( temp );
 
            output << "\n";
            output.close();
        }
    }

    static void output_length_differential(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        auto compares = get_smp_compare( bed_samples );
    }

    static void output_arms_differential(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        auto compares = get_smp_compare( bed_samples );
    }

};

} // end of namespace algorithm
} // end of namespace ago
