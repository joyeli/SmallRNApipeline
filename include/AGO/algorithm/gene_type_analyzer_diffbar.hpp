#pragma once
#include <cmath>
#include <boost/math/special_functions/beta.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDiffBar
{
  public:

    using TargetScanType = std::map< std::string, std::map< std::string, std::set< std::string >>>;
    //                                mir_name                mir_seed               targets

    GeneTypeAnalyzerDiffBar()
    {}

    static std::map< std::string, std::set< std::string >> get_anno_mapping( auto& bed_samples, auto& filter_ppm )
    {
        std::map< std::string, std::set< std::string >> anno_mapping;

        std::string mirbase = "";
        std::string gencode = "";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                if( raw_bed.annotation_info_.empty() || raw_bed.annotation_info_[0].empty() ) continue;
                if( raw_bed.annotation_info_[0][0] != "miRNA" ) continue;

                for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); i++ )
                    for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                        if( raw_bed.annotation_info_[i][j] == "mirbase" )
                            mirbase = raw_bed.annotation_info_[i][ j+1 ];

                if( mirbase == "" ) continue;

                gencode = raw_bed.annotation_info_[0][1];
                anno_mapping[ mirbase ].emplace( gencode );

                mirbase = "";
            }
        }

        return anno_mapping;
    }

    static std::string seed_u2t( const std::string& seed )
    {
        std::string res;
        for( auto& s : seed )
            if( s == 'U' )
                res += 'T';
            else
                res += s;
        return res;
    }

    static TargetScanType get_targetscan( auto& bed_samples, auto& targetscan_path, auto& species_code, auto& pct_cutoff, auto& filter_ppm )
    {
        TargetScanType targetscan;
        std::map< std::string, std::set< std::string >> anno_mapping = get_anno_mapping( bed_samples, filter_ppm );

        std::string line;
        std::string seed;

        std::vector< std::string > split;
        std::ifstream infile( targetscan_path );

        while( std::getline( infile, line ))
        {
            if( line.substr( 0, 4 ) != "ENST" ) continue;
            boost::iter_split( split, line, boost::algorithm::first_finder( "\t" ));

            if( species_code != "" && split[13].substr( 0, 3 ) != species_code ) continue;
            if( pct_cutoff != 0 && split[16] == "NULL" ) continue;
            if( pct_cutoff != 0 && std::stod( split[16] ) < pct_cutoff ) continue;
            if( anno_mapping.find( split[13] ) == anno_mapping.end() ) continue;

            seed = seed_u2t( split[2] );

            for( auto& gencode : anno_mapping[ split[13] ])
                targetscan[ gencode ][ seed ].emplace( split[1] + "\t" + split[16] );
        }

        return targetscan;
    }

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
            + s1 + ":" + s2 + ":Utail\t"
            + s1 + ":" + s2 + ":Other\t"
            + s2 + ":" + s1 + ":GMPM\t"
            + s2 + ":" + s1 + ":GM\t"
            + s2 + ":" + s1 + ":PM\t"
            + s2 + ":" + s1 + ":Atail\t"
            + s2 + ":" + s1 + ":Ctail\t"
            + s2 + ":" + s1 + ":Gtail\t"
            + s2 + ":" + s1 + ":Utail\t"
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
                case 'G' : i = 2; break;
                case 'U' : i = 3; break;
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
        double fold_change1 = s1 / s2;
        double fold_change2 = s2 / s1;
        double count1 = std::log( s1 );
        double count2 = std::log( s2 );
        // double count1 = s1 / ( std::log( s1 > s2 ? s1 : s2 ) / std::log(2) );
        // double count2 = s2 / ( std::log( s1 > s2 ? s1 : s2 ) / std::log(2) );
        double p_value = 1 - ( fold_change1 < 1
                ? boost::math::ibeta(( 3.26 + count2 ), ( 1.63 + count1 ), ( fold_change2 / ( 1 + fold_change2 )))
                : boost::math::ibeta(( 3.26 + count1 ), ( 1.63 + count2 ), ( fold_change1 / ( 1 + fold_change1 )))
                );
        return std::make_pair( fold_change1, p_value );
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
        {
            std::stringstream ss;
            ss << diff.second;
            out += "\t" + std::to_string( diff.first ) + ":" + ss.str();
        }

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
            std::vector< std::vector< CountingTableType >>& anno_table_tail_isomir,
            TargetScanType& targetscan,
            const bool& is_isomir = true
            )
    {
        std::ofstream output;

        //                      s1           s2                    type
        std::map< std::pair< std::string, std::string >, std::map< char,
            std::map< std::string, std::pair< double, double >>>> diffs;
        //               miR                   fold   pvalue

        for( auto& compare : get_smp_compare( bed_samples ))
        {
            auto& s1 = compare[0];
            auto& s2 = compare[1];

            output.open( output_path + "LoadingDifferential_" + s1.second + "_" + s2.second + ( is_isomir ? "-isomiRs" : "" ) + ".text" );
            output << output_header( s1.second, s2.second );

            double gm1, gm2, pm1, pm2, at1, at2, ct1, ct2, gt1, gt2, tt1, tt2, ot1, ot2;
            std::vector< std::pair< double, double >> diff_vec;
            std::vector< std::tuple< double, double, std::string >> out_temp;
            //                       total  smallestPvalue

            std::vector< std::string > split;
            std::set< std::string > anno_idx;
            std::vector< std::vector< CountingTableType >> anno_table_tail;

            if( !is_isomir )
            {
                anno_table_tail = std::vector< std::vector< CountingTableType >>( anno_table_tail_isomir.size(), std::vector< CountingTableType >() );

                for( std::size_t smp = 0; smp < anno_table_tail_isomir.size(); smp++ )
                    anno_table_tail[ smp ] = std::vector< CountingTableType >( anno_table_tail_isomir[ smp ].size(), CountingTableType() );

                for( auto& anno : ano_len_idx.first )
                {
                    boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
                    anno_idx.emplace( split[0] );

                    for( std::size_t smp = 0; smp < anno_table_tail_isomir.size(); smp++ )
                    {
                        for( std::size_t i = 0; i < 6; i++ )
                        {
                            if( anno_table_tail[ smp ][i].find( split[0] ) == anno_table_tail[ smp ][i].end() )
                                for( auto& len : ano_len_idx.second ) anno_table_tail[ smp ][i][ split[0] ][ len ] = 0.0;

                            for( auto& len : anno_table_tail_isomir[ smp ][i][ anno ] )
                                anno_table_tail[ smp ][i][ split[0] ][ len.first ] += anno_table_tail_isomir[ smp ][i][ anno ][ len.first ];
                        }
                    }
                }
            }
            else anno_idx = ano_len_idx.first;

            for( auto& anno : anno_idx )
            {
                gm1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'M' );
                gm2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'M' );
                pm1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'P' );
                pm2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'P' );
                at1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'A' );
                at2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'A' );
                ct1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'C' );
                ct2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'C' );
                gt1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'G' );
                gt2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'G' );
                tt1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'U' );
                tt2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'U' );
                ot1 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s1.first, anno, ano_len_idx.second, 'O' );
                ot2 = get_value(( is_isomir ? anno_table_tail_isomir : anno_table_tail ), s2.first, anno, ano_len_idx.second, 'O' );

                diff_vec.emplace_back( get_differential( gm1 + pm1 + 1, gm2 + pm2  + 1 ));
                diff_vec.emplace_back( get_differential( gm1       + 1, gm2        + 1 ));
                diff_vec.emplace_back( get_differential(       pm1 + 1,       pm2  + 1 ));
                diff_vec.emplace_back( get_differential(    at1    + 1,    at2     + 1 ));
                diff_vec.emplace_back( get_differential(    ct1    + 1,    ct2     + 1 ));
                diff_vec.emplace_back( get_differential(    gt1    + 1,    gt2     + 1 ));
                diff_vec.emplace_back( get_differential(    tt1    + 1,    tt2     + 1 ));
                diff_vec.emplace_back( get_differential(    ot1    + 1,    ot2     + 1 ));

                diffs[{ s1.second, s2.second }][ 'S' ][ anno ] = diff_vec[ diff_vec.size() - 8 ];
                diffs[{ s1.second, s2.second }][ 'M' ][ anno ] = diff_vec[ diff_vec.size() - 7 ];
                diffs[{ s1.second, s2.second }][ 'P' ][ anno ] = diff_vec[ diff_vec.size() - 6 ];
                diffs[{ s1.second, s2.second }][ 'A' ][ anno ] = diff_vec[ diff_vec.size() - 5 ];
                diffs[{ s1.second, s2.second }][ 'C' ][ anno ] = diff_vec[ diff_vec.size() - 4 ];
                diffs[{ s1.second, s2.second }][ 'G' ][ anno ] = diff_vec[ diff_vec.size() - 3 ];
                diffs[{ s1.second, s2.second }][ 'U' ][ anno ] = diff_vec[ diff_vec.size() - 2 ];

                diff_vec.emplace_back( get_differential( gm2 + pm2 + 1, gm1 + pm1  + 1 ));
                diff_vec.emplace_back( get_differential( gm2       + 1, gm1        + 1 ));
                diff_vec.emplace_back( get_differential(       pm2 + 1,       pm1  + 1 ));
                diff_vec.emplace_back( get_differential(    at2    + 1,    at1     + 1 ));
                diff_vec.emplace_back( get_differential(    ct2    + 1,    ct1     + 1 ));
                diff_vec.emplace_back( get_differential(    gt2    + 1,    gt1     + 1 ));
                diff_vec.emplace_back( get_differential(    tt2    + 1,    tt1     + 1 ));
                diff_vec.emplace_back( get_differential(    ot2    + 1,    ot1     + 1 ));

                diffs[{ s2.second, s1.second }][ 'S' ][ anno ] = diff_vec[ diff_vec.size() - 8 ];
                diffs[{ s2.second, s1.second }][ 'M' ][ anno ] = diff_vec[ diff_vec.size() - 7 ];
                diffs[{ s2.second, s1.second }][ 'P' ][ anno ] = diff_vec[ diff_vec.size() - 6 ];
                diffs[{ s2.second, s1.second }][ 'A' ][ anno ] = diff_vec[ diff_vec.size() - 5 ];
                diffs[{ s2.second, s1.second }][ 'C' ][ anno ] = diff_vec[ diff_vec.size() - 4 ];
                diffs[{ s2.second, s1.second }][ 'G' ][ anno ] = diff_vec[ diff_vec.size() - 3 ];
                diffs[{ s2.second, s1.second }][ 'U' ][ anno ] = diff_vec[ diff_vec.size() - 2 ];

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

        output_preference( output_path, is_isomir, targetscan, diffs );
    }

    static std::set< std::string > get_target( auto& targetscan, auto& mir )
    {
        std::vector< std::string > split;
        boost::iter_split( split, mir, boost::algorithm::first_finder( "|" ));
        boost::iter_split( split, split[0], boost::algorithm::first_finder( "_" ));
        return targetscan[ split[0] ][ split[1] ];
    }

    static void output_preference(
            const std::string& output_path,
            const bool& is_isomir,
            auto& targetscan,
            auto& diffs
            )
    {
        std::vector< std::pair< char, std::string >> types(
        {
              { 'S', "GMPM/" }
            , { 'M', "GM/" }
            , { 'P', "PM/" }
            , { 'A', "Atail/" }
            , { 'C', "Ctail/" }
            , { 'G', "Gtail/" }
            , { 'U', "Utail/" }
        });

        boost::filesystem::create_directory( boost::filesystem::path( output_path + "../Preference/" ));

        for( auto& type : types )
            boost::filesystem::create_directory( boost::filesystem::path( output_path + "../Preference/" + type.second ));

        std::size_t smpidx = 0;
        std::set< std::string > samples;

        std::map< std::string, std::map< std::string, std::set< std::string >>> mir_unilike;
        std::map< std::string, std::map< std::string, std::set< std::string >>> mir_dislike;
        //          type                    mir                     samples

        std::map< std::string, std::map< std::string, std::set< std::string >>> set_unitgsc;
        std::map< std::string, std::map< std::string, std::set< std::string >>> set_distgsc;
        
        std::map< std::string, std::map< std::string, std::ofstream >> out_unilike;
        std::map< std::string, std::map< std::string, std::ofstream >> out_dislike;
        //          type                    sample      ofstream
        std::map< std::string, std::map< std::string, std::ofstream >> out_unitgsc;
        std::map< std::string, std::map< std::string, std::ofstream >> out_distgsc;

        for( auto& smp : diffs )
        {
            samples.emplace( smp.first.first );
            samples.emplace( smp.first.second );
        }

        for( auto& smp : samples ) for( auto& type : types )
        {
            out_unilike[ type.second ][ smp ].open( output_path + "../Preference/" + type.second + smp + "_UniLike" + ( is_isomir ? "-isomiRs" : "" ) + ".text" );
            out_dislike[ type.second ][ smp ].open( output_path + "../Preference/" + type.second + smp + "_DisLike" + ( is_isomir ? "-isomiRs" : "" ) + ".text" );

            if( !targetscan.empty() && is_isomir )
            {
                out_unitgsc[ type.second ][ smp ].open( output_path + "../Preference/" + type.second + smp + "_UniLike-isomiRs_target.text" );
                out_distgsc[ type.second ][ smp ].open( output_path + "../Preference/" + type.second + smp + "_DisLike-isomiRs_target.text" );
            }

            for( auto& compare : diffs )
            {
                smpidx = 0;

                if( compare.first.first  == smp ) smpidx = 1;
                if( compare.first.second == smp ) smpidx = 2;

                if( smpidx == 0 ) continue;

                for( auto& mir : compare.second[ type.first ])
                {
                    if( mir.second.second > 0.05 ) continue;
                    if( mir.second.first < 2 && mir.second.first > 0.5 ) continue;

                    if( smpidx == 1 && mir.second.first >= 2   ) mir_unilike[ type.second ][ mir.first ].emplace( compare.first.second );
                    if( smpidx == 1 && mir.second.first <= 0.5 ) mir_dislike[ type.second ][ mir.first ].emplace( compare.first.second );
                    if( smpidx == 2 && mir.second.first >= 2   ) mir_dislike[ type.second ][ mir.first ].emplace( compare.first.first  );
                    if( smpidx == 2 && mir.second.first <= 0.5 ) mir_unilike[ type.second ][ mir.first ].emplace( compare.first.first  );
                }
            }
        }

        for( auto& type : mir_unilike )
            for( auto& mir : type.second )
                if( mir.second.size() == samples.size() -1 )
                    for( auto& smp : samples )
                        if( mir.second.find( smp ) == mir.second.end() )
                        {
                            out_unilike[ type.first ][ smp ] << mir.first << "\n";

                            if( !targetscan.empty() && is_isomir )
                                for( auto& target : get_target( targetscan, mir.first ))
                                    set_unitgsc[ type.first ][ smp ].emplace( target );
                        }

        for( auto& type : mir_dislike )
            for( auto& mir : type.second )
                if( mir.second.size() == samples.size() -1 )
                    for( auto& smp : samples )
                        if( mir.second.find( smp ) == mir.second.end() )
                        {
                            out_dislike[ type.first ][ smp ] << mir.first << "\n";

                            if( !targetscan.empty() && is_isomir )
                                for( auto& target : get_target( targetscan, mir.first ))
                                    set_distgsc[ type.first ][ smp ].emplace( target );
                        }

        if( !targetscan.empty() && is_isomir )
        {
            for( auto& type : set_unitgsc )
                for( auto& smp : type.second )
                    for( auto& target : smp.second )
                        out_unitgsc[ type.first ][ smp.first ] << target << "\n";

            for( auto& type : set_distgsc )
                for( auto& smp : type.second )
                    for( auto& target : smp.second )
                        out_distgsc[ type.first ][ smp.first ] << target << "\n";
        }

        for( auto& smp : samples ) for( auto& type : types )
        {
            out_unilike[ type.second ][ smp ].close();
            out_dislike[ type.second ][ smp ].close();

            if( !targetscan.empty() && is_isomir )
            {
                out_unitgsc[ type.second ][ smp ].close();
                out_distgsc[ type.second ][ smp ].close();
            }
        }
    }

    static void output_diffbar_isomirs(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& biotype,
            const std::string& token
            )
    {
        std::ofstream output;
        std::vector< std::pair< std::string, double >> total_vec;
        std::map< std::string, std::map< std::string, std::vector< double >>> mirlen_map;
        GeneTypeAnalyzerBarplot::make_table_isomirs( bed_samples, ano_len_idx, anno_table_tail, token, total_vec, mirlen_map );
        GeneTypeAnalyzerBarplot::sort_total_vec( total_vec );

        for( auto& len : ano_len_idx.second )
        {
            output.open( output_name + token + "_" + std::to_string( len ) + "-isomiRs.tsv" );
            GeneTypeAnalyzerBarplot::output_table( output, std::to_string( len ), bed_samples, total_vec, mirlen_map );
            output.close();
        }
    }

    static void output_diffbar(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& biotype,
            const std::string& token
            )
    {
        std::ofstream output;
        std::vector< std::pair< std::string, double >> total_vec;
        std::map< std::string, std::map< std::string, std::vector< double >>> mirlen_map;
        GeneTypeAnalyzerBarplot::make_table( bed_samples, ano_len_idx, anno_table_tail, token, total_vec, mirlen_map );
        GeneTypeAnalyzerBarplot::sort_total_vec( total_vec );

        for( auto& len : ano_len_idx.second )
        {
            output.open( output_name + token + "_" + std::to_string( len ) + ".tsv" );
            GeneTypeAnalyzerBarplot::output_table( output, std::to_string( len ), bed_samples, total_vec, mirlen_map );
            output.close();
        }
    }

    static void output_diffbar_visualization( const std::string& output_name, const bool& isSeed )
    {
        if( !boost::filesystem::exists( output_name + "GMPM.tsv" ))
             boost::filesystem::create_symlink( "../ValPlot/GMPM.tsv", ( output_name + "GMPM.tsv" ).c_str() );

        if( !isSeed && !boost::filesystem::exists( output_name + "GMPM-isomiRs.tsv" ))
             boost::filesystem::create_symlink( "../ValPlot/GMPM-isomiRs.tsv", ( output_name + "GMPM-isomiRs.tsv" ).c_str() );

        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $Diff = $_POST['Diff'];" << "\n";
        output << "        $GMPM = $_POST['GMPM'];" << "\n";
        output << "        $is5p3p = $_POST['is5p3p'];" << "\n";
        output << "        $IsomiRs = $_POST['IsomiRs'];" << "\n";
        output << "        $P_Value = $_POST['P_Value'];" << "\n";
        output << "        $Sample1 = $_POST['Sample1'];" << "\n";
        output << "        $Sample2 = $_POST['Sample2'];" << "\n";
        output << "        $MaxHight = $_POST['MaxHight'];" << "\n";
        output << "        $Top_miRNA = $_POST['Top_miRNA'];" << "\n";
        output << "        $Min_Length = $_POST['Min_Length'];" << "\n";
        output << "        $Max_Length = $_POST['Max_Length'];" << "\n";
        output << "        $Fold_Change = $_POST['Fold_Change'];" << "\n";
        output << "        $barPlotType = $_POST['barPlotType'];" << "\n";
        output << "        $Select_Type = $_POST['Select_Type'];" << "\n";
        output << "        $Single_Anno = $_POST['Single_Anno'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js></script>';" << "\n";
        output << "        echo '<link href=http://192.168.1.11:6680/WorkDir/AgoD3/ForAgoSorting/lib/svg0331.css rel=stylesheet type=text/css>';" << "\n";
        output << "        echo '<style> .x.axis path { display: none; }</style>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Diff ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Diff onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '.( $Diff=='' ? 'selected' : '' ).'>Diff By</option>';" << "\n";
        output << "" << "\n";
        output << "        $Diff_List = array('GMPM', 'GM', 'PM', 'Atail', 'Ctail', 'Gtail', 'Utail', 'Other');" << "\n";
        output << "" << "\n";
        output << "        $Diff_Size = Count( $Diff_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Diff_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Diff_List[$i].' ';" << "\n";
        output << "            echo $Diff == $Diff_List[$i] ? 'selected ' : '';" << "\n";
        output << "            echo '>' . $Diff_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=GMPM onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '.( $GMPM=='' ? 'selected' : '' ).'>GM or PM</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPM_List = array('GMPM', 'GM', 'PM', 'Tailing', 'Atail', 'Ctail', 'Gtail', 'Utail', 'Other');" << "\n";
        output << "" << "\n";
        output << "        $GMPM_Size = Count( $GMPM_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $GMPM_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$GMPM_List[$i].' ';" << "\n";
        output << "            echo $GMPM == $GMPM_List[$i] ? 'selected ' : '';" << "\n";
        output << "            echo '>' . $GMPM_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .text' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        $Sample_List = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Temp = Explode( '_', $TSV_List[$i] );" << "\n";
        output << "            Array_Push( $Sample_List, $Temp[1] );" << "\n";
        output << "" << "\n";
        output << "            $Temp = Explode( '.', $Temp[2] );" << "\n";
        output << "            $Temp = Explode( '-', $Temp[0] );" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $Sample_List, $Temp[0] );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Sample_List = Array_Unique( $Sample_List );" << "\n";
        output << "        $Sample_List = Array_Values( $Sample_List );" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Sample1 =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Sample1 onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Sample1=='') echo 'selected'; echo '>Select Sample1</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Sample_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Sample_List[$i] == '' ) continue;" << "\n";
        output << "            if( $Sample_List[$i] == $Sample2 ) continue;" << "\n";
        output << "            echo '<option value='.$Sample_List[$i].' ';" << "\n";
        output << "            if( $Sample1 == $Sample_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Sample_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Sample2 =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Sample2 onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Sample2=='') echo 'selected'; echo '>Select Sample2</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Sample_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Sample_List[$i] == '' ) continue;" << "\n";
        output << "            if( $Sample_List[$i] == $Sample1 ) continue;" << "\n";
        output << "            echo '<option value='.$Sample_List[$i].' ';" << "\n";
        output << "            if( $Sample2 == $Sample_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Sample_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";

        if( !isSeed )
        {
            output << "#<!--================== IsomiRs =====================-->" << "\n";
            output << "                " << "\n";
            output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
            output << "" << "\n";
            output << "        echo '<select name=IsomiRs onchange=this.form.submit();>';" << "\n";
            output << "        echo '<option '; if($IsomiRs=='') echo 'selected'; echo '>Show IsomiRs?</option>';" << "\n";
            output << "" << "\n";
            output << "        $miR_List = array('Yes', 'No');" << "\n";
            output << "" << "\n";
            output << "        For( $i = 0; $i < Count( $miR_List ); ++$i )" << "\n";
            output << "        {" << "\n";
            output << "            echo '<option value='.$miR_List[$i].' ';" << "\n";
            output << "" << "\n";
            output << "            if( $IsomiRs == $miR_List[$i] )" << "\n";
            output << "                echo 'selected ';" << "\n";
            output << "" << "\n";
            output << "            echo '>' . $miR_List[$i] . '</option>';" << "\n";
            output << "        }" << "\n";
            output << "" << "\n";
            output << "        echo \"</select>" << "\n";
            output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
            output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
            output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
            output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
            output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
            output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
            output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
            output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
            output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
            output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
            output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
            output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
            output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
            output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
            output << "            </form>\";" << "\n";
        }
        else
        {
            output << "        $IsomiRs = 'No';" << "\n";
        }

        output << "" << "\n";
        output << "#<!--================ 5p3pSelector ==================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=is5p3p onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($is5p3p=='') echo 'selected'; echo 'value= >5p3p</option>';" << "\n";
        output << "" << "\n";
        output << "        $Arm_List = array( '5p', '3p' );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Arm_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Arm_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $is5p3p == $Arm_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Arm_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== barPlotType ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=barPlotType onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $barPlotType_List = array('ppm', 'log2', '100%');" << "\n";
        output << "        $barPlotType_Size = Count( $barPlotType_List );" << "\n";
        output << "        if( $barPlotType == '' ) $barPlotType = $barPlotType_List[0];" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $barPlotType_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$barPlotType_List[$i].' ';" << "\n";
        output << "            echo $barPlotType == $barPlotType_List[$i] ? 'selected ' : '';" << "\n";
        output << "            echo '>' . $barPlotType_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Select Type ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Select_Type onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $Select_Type_List = array('SelectType', 'SingleAnno', 'TopRanking');" << "\n";
        output << "        $Select_Type_Size = Count( $Select_Type_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Select_Type_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Select_Type_List[$i].' ';" << "\n";
        output << "            echo $Select_Type == $Select_Type_List[$i] ? 'selected ' : '';" << "\n";
        output << "            echo '>' . $Select_Type_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Top miRNA Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Select_Type == 'TopRanking' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=Top_miRNA size=4 value=';" << "\n";
        output << "            echo $Top_miRNA == '' ? 'miRNA#' : $Top_miRNA;" << "\n";
        output << "            echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "                <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "                <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "                <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "                <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "                <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "                <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "                <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "                <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "                <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "                <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "                <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "                <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "                <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "                <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Read File ====================-->" << "\n";
        output << "        " << "\n";
        output << "        $TSV_File = '';" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $Total_PPM = 0;" << "\n";
        output << "" << "\n";
        output << "        if( File_Exists( './LoadingDifferential_'.$Sample1.'_'.$Sample2.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text' ))" << "\n";
        output << "            $TSV_File =  './LoadingDifferential_'.$Sample1.'_'.$Sample2.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text';" << "\n";
        output << "" << "\n";
        output << "        if( File_Exists( './LoadingDifferential_'.$Sample2.'_'.$Sample1.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text' ))" << "\n";
        output << "            $TSV_File =  './LoadingDifferential_'.$Sample2.'_'.$Sample1.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text';" << "\n";
        output << "" << "\n";
        output << "        $Filtered_miRNAs = Array();" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                for( $i = 2; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    if( $inFile_Line[$i] == $Sample1.':'.$Sample2.':'.$Diff )" << "\n";
        output << "                        $Column = $i;" << "\n";
        output << "" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Value = Explode( ':', $inFile_Line[ $Column ] );" << "\n";
        output << "" << "\n";
        output << "            if( $Fold_Change != '' && $Fold_Change != '> |Log2( Fold )|' && Abs( Log( $Value[0], 2 )) < $Fold_Change ) continue;" << "\n";
        output << "            if( $P_Value     != '' && $P_Value     != '< P-Value'        &&           $Value[1]       > $P_Value     ) continue;" << "\n";
        output << "" << "\n";
        output << "            if( $Total_PPM < $inFile_Line[1] ) $Total_PPM = $inFile_Line[1];" << "\n";
        output << "" << "\n";
        output << "            $Filtered_miRNAs[ $inFile_Line[0] ] =" << "\n";
        output << "                [" << "\n";
        output << "                    Number_format( (float)$inFile_Line[1]     , 2, '.', '' )," << "\n";
        output << "                    Number_format( (float)Log(  $Value[0], 2 ), 2, '.', '' )," << "\n";
        output << "                    Number_format( (float)      $Value[1]     , 3, '.', '' )" << "\n";
        output << "                ];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $Total_PPM = 100 / $Total_PPM;" << "\n";
        output << "" << "\n";
        output << "        $Single_Anno_List = Array();" << "\n";
        output << "        $inFile = new SplFileObject( './GMPM'.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.tsv' );" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "            $inFile_Line = Explode( '!', $inFile_Line[0] );" << "\n";
        output << "            $inFile_Line = Explode( '*', $inFile_Line[0] );" << "\n";
        output << "" << "\n";
        output << "            if( !Array_Key_Exists( $inFile_Line[0], $Filtered_miRNAs ))" << "\n";
        output << "                continue;" << "\n";
        output << "" << "\n";
        output << "            if( $is5p3p == '5p' && StrPos( $inFile_Line[0], '5p' ) === false ) continue;" << "\n";
        output << "            if( $is5p3p == '3p' && StrPos( $inFile_Line[0], '3p' ) === false ) continue;" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $Single_Anno_List, $inFile_Line[0] );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = $Min_Length; $i <= $Max_Length; ++$i )" << "\n";
        output << "            Array_Push( $TSV_File, './'.$GMPM.'_'.$i.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.tsv' );" << "\n";
        output << "" << "\n";
        output << "#<!--================== Single Annotation ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Select_Type == 'SingleAnno' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<select name=Single_Anno onchange=this.form.submit();>';" << "\n";
        output << "            echo '<option '.( $Single_Anno == '' ? 'selected' : '' ).'>Annotation (ppm, log2Fold, Pvalue)</option>';" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Single_Anno_List as $idx => $anno )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$anno.' ';" << "\n";
        output << "                echo $Single_Anno == $anno ? 'selected ' : '';" << "\n";
        output << "                echo '>'.$anno.' ('.$Filtered_miRNAs[ $anno ][0].', '.$Filtered_miRNAs[ $anno ][1].', '.$Filtered_miRNAs[ $anno ][2].')</option>';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>" << "\n";
        output << "                <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "                <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "                <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "                <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "                <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "                <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "                <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "                <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "                <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "                <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "                <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "                <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "                <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>var select_color_map = d3.scale.linear().domain([ 0, 100 ]).range([ 'WhiteSmoke', 'Black' ]);\";" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Single_Anno_List as $idx => $anno )" << "\n";
        output << "            {" << "\n";
        output << "               echo \"$( 'option[value=\\\"$anno\\\"]' ).css(" << "\n";
        output << "               {" << "\n";
        output << "                   'background-color': select_color_map( '\".( $Filtered_miRNAs[ $anno ][0] * $Total_PPM ).\"' )," << "\n";
        output << "                   'color': '\".(( $Filtered_miRNAs[ $anno ][0] * $Total_PPM ) > 25 ? 'white' : 'black'  ).\"'" << "\n";
        output << "               });\";" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '</script>';" << "\n";
        output << "        }" << "\n";
        output << "        else $Single_Anno = '';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Max Hight ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $barPlotType != '100%' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<input type=text name=MaxHight size=8 value=';" << "\n";
        output << "            echo $MaxHight == '' ? 'MaxHight' : $MaxHight;" << "\n";
        output << "            echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MinMax Length =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Min_Length size=3 value=';" << "\n";
        output << "        echo $Min_Length == '' ? 'minLen' : $Min_Length;" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Max_Length size=3 value=';" << "\n";
        output << "        echo $Max_Length == '' ? 'maxLen' : $Max_Length;" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Filter Fold Pvalue =====================-->" << "\n";
        output << "        " << "\n";
        output << "        echo '<input type=text name=Fold_Change size=10 value=';" << "\n";
        output << "        echo $Fold_Change == '' || $Fold_Change == '> |Log2( Fold )|'? '\"> |Log2( Fold )|\"' : $Fold_Change;" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=P_Value size=5 value=';" << "\n";
        output << "        echo $P_Value == '' || $P_Value == '< P-Value' ? '\"< P-Value\"' : $P_Value;" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Switch =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='Diff' value='$Diff' />" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='is5p3p' value='$is5p3p' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='P_Value' value='$P_Value' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='Top_miRNA' value='$Top_miRNA' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='Fold_Change' value='$Fold_Change' />" << "\n";
        output << "            <input type='hidden' name='barPlotType' value='$barPlotType' />" << "\n";
        output << "            <input type='hidden' name='Select_Type' value='$Select_Type' />" << "\n";
        output << "            <input type='hidden' name='Single_Anno' value='$Single_Anno' />" << "\n";
        output << "            <input type='submit' value='Switch Samples' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Read Multi TSV ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( Count( $Filtered_miRNAs ) == 0 )" << "\n";
        output << "            echo '<h3>No Data</h3><br/>';" << "\n";
        output << "" << "\n";
        output << "        $MaxValue = 0;" << "\n";
        output << "" << "\n";
        output << "        Foreach( $TSV_File as $idx => $tsv )" << "\n";
        output << "        {" << "\n";
        output << "            if( !File_Exists( $tsv )) continue;" << "\n";
        output << "" << "\n";
        output << "            $isHeader = true;" << "\n";
        output << "            $inFile = new SplFileObject( $tsv );" << "\n";
        output << "" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "                if( $isHeader )" << "\n";
        output << "                {" << "\n";
        output << "                    $isHeader = false;" << "\n";
        output << "                    continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "                if( $Single_Anno != '' && $inFile_Line[0] != $Single_Anno ) continue;" << "\n";
        output << "                if( !Array_Key_Exists( $inFile_Line[0], $Filtered_miRNAs )) continue;" << "\n";
        output << "" << "\n";
        output << "                For( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    if( $inFile_Line[$i] >= $MaxValue ) $MaxValue = $inFile_Line[$i];" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Temp_Array = Array();" << "\n";
        output << "        $MaxValue = $barPlotType != 'ppm' && $barPlotType != 'log2' ? 100 : $MaxValue;" << "\n";
        output << "" << "\n";
        output << "        $MaxValue = $barPlotType == 'ppm'  ?    ( $MaxHight == '' || $MaxHight == 'MaxHight' ? $MaxValue : $MaxHight    ) : $MaxValue;" << "\n";
        output << "        $MaxValue = $barPlotType == 'log2' ? Log( $MaxHight == '' || $MaxHight == 'MaxHight' ? $MaxValue : $MaxHight, 2 ) : $MaxValue;" << "\n";
        output << "" << "\n";
        output << "        Foreach( $TSV_File as $idx => $tsv )" << "\n";
        output << "        {" << "\n";
        output << "            if( !File_Exists( $tsv )) continue;" << "\n";
        output << "" << "\n";
        output << "            $count = 0;" << "\n";
        output << "            $isHeader = true;" << "\n";
        output << "            $inFile = new SplFileObject( $tsv );" << "\n";
        output << "" << "\n";
        output << "            $Temp_Array[ $idx ] = Tempnam( '/tmp', $tsv );" << "\n";
        output << "            $Ftemp = fopen( $Temp_Array[ $idx ], 'w' );" << "\n";
        output << "" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "                if( $isHeader )" << "\n";
        output << "                {" << "\n";
        output << "                    $isHeader = false;" << "\n";
        output << "                    fwrite( $Ftemp, $inFile_Lines );" << "\n";
        output << "                    continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "                if( $is5p3p == '5p' && StrPos( $inFile_Line[0], '5p' ) === false ) continue;" << "\n";
        output << "                if( $is5p3p == '3p' && StrPos( $inFile_Line[0], '3p' ) === false ) continue;" << "\n";
        output << "" << "\n";
        output << "                if( $Single_Anno != '' && $inFile_Line[0] != $Single_Anno ) continue;" << "\n";
        output << "                if( !Array_Key_Exists( $inFile_Line[0], $Filtered_miRNAs )) continue;" << "\n";
        output << "" << "\n";
        output << "                if( $barPlotType != 'ppm' && $barPlotType != 'log2' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Total_Value = 0;" << "\n";
        output << "                    fwrite( $Ftemp, $inFile_Line[0] );" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 1; $i < Count( $inFile_Line ); ++$i ) $Total_Value += $inFile_Line[ $i ];" << "\n";
        output << "                    For( $i = 1; $i < Count( $inFile_Line ); ++$i ) fwrite( $Ftemp, \"\\t\".Number_format( $inFile_Line[ $i ] * 100 / $Total_Value, 0 ));" << "\n";
        output << "" << "\n";
        output << "                    fwrite( $Ftemp, \"\\n\");" << "\n";
        output << "                }" << "\n";
        output << "                else if( $barPlotType == 'log2' )" << "\n";
        output << "                {" << "\n";
        output << "                    fwrite( $Ftemp, $inFile_Line[0] );" << "\n";
        output << "                    For( $i = 1; $i < Count( $inFile_Line ); ++$i ) fwrite( $Ftemp, \"\\t\".Number_format( $inFile_Line[ $i ] == 0 ? 0 : Log( $inFile_Line[ $i ], 2 )));" << "\n";
        output << "                    fwrite( $Ftemp, \"\\n\");" << "\n";
        output << "                }" << "\n";
        output << "                else fwrite( $Ftemp, $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "                $count++;" << "\n";
        output << "                if( $count >= $Top_miRNA ) break;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            fclose( $Ftemp );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== barPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        //========== svg view var set ==========" << "\n";
        output << "" << "\n";
        output << "        $width = ( $Top_miRNA == '' ? 1 : $Top_miRNA ) * 180;" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "                var margin = {top: 20, right: 60, bottom: 35, left: 40}," << "\n";
        output << "                    width = $width - margin.left - margin.right, //180(per one miRNA)" << "\n";
        output << "                    height = 300 - margin.top - margin.bottom;" << "\n";
        output << "                " << "\n";
        output << "                var x0 = d3.scale.ordinal()" << "\n";
        output << "                    .rangeRoundBands([0, width], .3);" << "\n";
        output << "                " << "\n";
        output << "                var x1 = d3.scale.ordinal();" << "\n";
        output << "                " << "\n";
        output << "                var y = d3.scale.linear()" << "\n";
        output << "                    .range([height, 0]);" << "\n";
        output << "                " << "\n";
        output << "                var color = d3.scale.ordinal()" << "\n";
        output << "                    .range(['#F75C2F', '#E8B647', '#838A2D', '#66BAB7', '#6E75A4', '#72636E']);" << "\n";
        output << "                " << "\n";
        output << "                var xAxis = d3.svg.axis()" << "\n";
        output << "                    .scale(x0)" << "\n";
        output << "                    .orient('bottom');" << "\n";
        output << "                " << "\n";
        output << "                var yAxis = d3.svg.axis()" << "\n";
        output << "                    .scale(y)" << "\n";
        output << "                    .orient('left')" << "\n";
        output << "                    .tickFormat(d3.format('.2s'));\";" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Temp_Array as $idx => $tsv )" << "\n";
        output << "        {" << "\n";
        output << "            echo \"var svg$idx = d3.select('body').append('svg')" << "\n";
        output << "                    .attr('id', 'svg$idx')" << "\n";
        output << "                .attr('width', width + margin.left + margin.right)" << "\n";
        output << "                .attr('height', height + margin.top + margin.bottom)" << "\n";
        output << "                .append('g')" << "\n";
        output << "                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n";
        output << "            " << "\n";
        output << "                d3.tsv('$tsv', function(error, data) {" << "\n";
        output << "                var sample = d3.keys(data[0]).filter(function(key) { return key !== 'Annotation'; });" << "\n";
        output << "            " << "\n";
        output << "                data.forEach(function(d) {" << "\n";
        output << "                     d.values = sample.map(function(name) { return {name: name, value: +d[name]}; });" << "\n";
        output << "                });" << "\n";
        output << "    " << "\n";
        output << "                x0.domain(data.map(function(d) { return d.Annotation; }));" << "\n";
        output << "                x1.domain(sample).rangeRoundBands([0, x0.rangeBand()]);" << "\n";
        output << "                y.domain([0, $MaxValue]);" << "\n";
        output << "            " << "\n";
        output << "                    svg$idx.append('g')" << "\n";
        output << "                    .attr('class', 'x axis')" << "\n";
        output << "                    .attr('transform', 'translate(0,' + height + ')')" << "\n";
        output << "                    .call(xAxis)" << "\n";
        output << "                    .selectAll('text')" << "\n";
        output << "                    .style('font-size', function(d){" << "\n";
        output << "                        return 75 - ( d.length > 15 ? (( d.length - 15 ) * 1.25 ) : 0 ) + '%';" << "\n";
        output << "                    });" << "\n";
        output << "            " << "\n";
        output << "                    svg$idx.append('g')" << "\n";
        output << "                    .attr('class', 'y axis')" << "\n";
        output << "                    .call(yAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('transform', 'rotate(-90)')" << "\n";
        output << "                    .attr('y', 6)" << "\n";
        output << "                    .attr('dy', '.71em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                        .text('$TSV_File[$idx]');" << "\n";
        output << "            " << "\n";
        output << "                    var barchat = svg$idx.selectAll('.barchat')" << "\n";
        output << "                    .data(data).enter()" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', 18)" << "\n";
        output << "                    .attr('y1', 0)" << "\n";
        output << "                    .attr('x2', 18)" << "\n";
        output << "                    .attr('y2', height)" << "\n";
        output << "                    .attr('style', 'stroke:rgb(192,192,192);stroke-width:1')" << "\n";
        output << "                    .attr('transform', function(d) { return 'translate(' + (x0(d.Annotation) + d.values.length * x1.rangeBand()) + ',0)'; });" << "\n";
        output << "            " << "\n";
        output << "                    barchat = svg$idx.selectAll('.barchat')" << "\n";
        output << "                    .data(data)" << "\n";
        output << "                    .enter().append('g')" << "\n";
        output << "                    .attr('class', 'g')" << "\n";
        output << "                    .attr('transform', function(d) { return 'translate(' + x0(d.Annotation) + ',0)'; });" << "\n";
        output << "    " << "\n";
        output << "                var div = d3.select('body').append('div')" << "\n";
        output << "                    .attr('class', 'tooltip')" << "\n";
        output << "                    .attr('id', 'Mir')" << "\n";
        output << "                    .style('opacity', 0);" << "\n";
        output << "            " << "\n";
        output << "                barchat.selectAll('rect')" << "\n";
        output << "                    .data(function(d) { return d.values; })" << "\n";
        output << "                    .enter().append('rect')" << "\n";
        output << "                    .attr('width', x1.rangeBand() - 5)" << "\n";
        output << "                    .attr('x', function(d) { return x1(d.name); })" << "\n";
        output << "                    .attr('y', function(d){" << "\n";
        output << "                        if( d.value < $MaxValue )" << "\n";
        output << "                            return y(d.value);" << "\n";
        output << "                        else" << "\n";
        output << "                            return y($MaxValue);" << "\n";
        output << "                    })" << "\n";
        output << "                    .attr('height', function(d){" << "\n";
        output << "                        if( d.value < $MaxValue )" << "\n";
        output << "                            return height - y(d.value);" << "\n";
        output << "                        else" << "\n";
        output << "                            return height - y($MaxValue);" << "\n";
        output << "                    })" << "\n";
        output << "                    .style('fill', function(d) { return color(d.name); })" << "\n";
        output << "                    .on('mouseover', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(200)" << "\n";
        output << "                            .style('opacity', .9);" << "\n";
        output << "    " << "\n";
        output << "                        div.html('<table align=center ><tr><th>Sample</th><th>Value</th></tr><tr><th>' +" << "\n";
        output << "                                 d.name + '</th><th>' + d.value + '</th></tr>')" << "\n";
        output << "                            .style('left', (d3.event.pageX) + 'px')" << "\n";
        output << "                            .style('top', (d3.event.pageY - 28) + 'px');" << "\n";
        output << "                    })" << "\n";
        output << "                    .on('mouseout', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(500)" << "\n";
        output << "                            .style('opacity', 0);" << "\n";
        output << "                    });" << "\n";
        output << "    " << "\n";
        output << "                    var legend = svg$idx.selectAll('.legend')" << "\n";
        output << "                    .data(sample.slice())" << "\n";
        output << "                    .enter().append('g')" << "\n";
        output << "                    .attr('class', 'legend')" << "\n";
        output << "                    .attr('transform', function(d, i) { return 'translate(0,' + i * 20 + ')'; });" << "\n";
        output << "            " << "\n";
        output << "                legend.append('rect')" << "\n";
        output << "                    .attr('x', width + 40)" << "\n";
        output << "                    .attr('width', 18)" << "\n";
        output << "                    .attr('height', 18)" << "\n";
        output << "                    .style('fill', color);" << "\n";
        output << "            " << "\n";
        output << "                legend.append('text')" << "\n";
        output << "                    .attr('x', width + 34)" << "\n";
        output << "                    .attr('y', 9)" << "\n";
        output << "                    .attr('dy', '.35em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text(function(d) { return d; });" << "\n";
        output << "            });\";" << "\n";
        output << "        }" << "\n";
        output << "        echo '</script>';" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    } 
};

} // end of namespace algorithm
} // end of namespace ago
