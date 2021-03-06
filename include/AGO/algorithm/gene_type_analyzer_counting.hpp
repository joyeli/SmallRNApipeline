#pragma once
#include <AGO/format/md_rawbed.hpp>
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerCounting
{
  public:

    GeneTypeAnalyzerCounting()
    {}

    static AnnoLengthIndexType get_ano_len_idx(
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            const std::size_t& filter_ppm,
            const std::string& biotype = "", // default annotation type
            const std::string& token = ""
            )
    {
        std::vector< std::set< std::string >> smp_ano_idx( bed_samples.size(), std::set< std::string >() );
        std::vector< std::set< std::size_t >> smp_len_idx( bed_samples.size(), std::set< std::size_t >() );
        std::vector< std::map< std::string, double >> anno_ppm_check( bed_samples.size(), std::map< std::string, double >() );

        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &genome_table, &biotype, &bed_samples, &smp_ano_idx, &anno_ppm_check, &smp_len_idx, &filter_ppm, token ] () mutable
            {
                std::string gene_name;
                std::string gene_seed;
                
                bool isbiotype = false;

                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    if( raw_bed.ppm_ < filter_ppm ) continue;
                    isbiotype = biotype == "" ? true : false;
                    smp_len_idx[ smp ].emplace( raw_bed.length_ - raw_bed.tail_length_ );

                    // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
                    {
                        std::size_t i = 0; // do first priority
                        if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                        {
                            if( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))
                                isbiotype = true;

                            if( biotype != "" && raw_bed.annotation_info_[i][0] == biotype )
                                isbiotype = true;

                            if( !isbiotype ) continue;

                            for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                            {
                                gene_seed = GeneTypeAnalyzerCounting::seqT2U( raw_bed.getReadSeq( genome_table ).substr( 1, 7 ))
                                        + ( GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag ) != "" ? ( "|" + GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag )) : "" );

                                gene_name = token == "Seed" ? gene_seed : ( biotype == ""
                                    ?   raw_bed.annotation_info_[i][j]
                                    : ( raw_bed.annotation_info_[i][ j+1 ] + "_" + gene_seed )
                                );

                                if(( token == "3p" || token == "5p" ) && get_arm( raw_bed.annotation_info_[i][ j+1 ] ) != token ) continue;

                                smp_ano_idx[ smp ].emplace( gene_name );

                                if( biotype != "" )
                                {
                                    if( anno_ppm_check[ smp ].find( gene_name ) == anno_ppm_check[ smp ].end() )
                                        anno_ppm_check[ smp ][ gene_name ] = 0.0;

                                    anno_ppm_check[ smp ][ gene_name ] += raw_bed.ppm_;
                                }
                            }
                        }
                    }
                }
            });
        }

        smp_parallel_pool.flush_pool();

        std::set< std::string > ano_idx;
        std::set< std::size_t > len_idx;

        std::set< std::string > ano_temp;
        std::set< std::size_t > len_temp;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            ano_temp.insert( smp_ano_idx[ smp ].begin(), smp_ano_idx[ smp ].end() );
            len_temp.insert( smp_len_idx[ smp ].begin(), smp_len_idx[ smp ].end() );
        }

        double sum = 0.0;

        if( biotype != "" )
        {
            for( auto& anno : ano_temp )
            {
                sum = 0.0;

                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                {
                    if( anno_ppm_check[ smp ].find( anno ) != anno_ppm_check[ smp ].end() )
                        sum += anno_ppm_check[ smp ][ anno ];
                }

                if( sum >= 1 ) ano_idx.emplace( anno );
            }
        }

        std::size_t count = 0;
        std::size_t last_len;

        for( auto& len : len_temp )
        {
            if( count == 0 ) last_len = len -1;
            if( len - last_len > 1 )
                for( std::size_t i = last_len + 1; i < len; ++i )
                    len_idx.emplace( i ); 
            count++;
        }

        return std::make_pair( ( biotype != "" ? ano_idx : ano_temp ), len_idx );
    }

    static std::size_t which_tail( const std::string& tail )
    {
        std::set< char > nt_check;
        for( auto& nt : tail ) nt_check.emplace( nt );

        if( nt_check.size() == 0 ) return 5; // GM
        if( nt_check.size() >  1 ) return 4; // Other_Tail

        switch( *nt_check.begin() )
        {
            case 'A' : return 0; // A_Tail
            case 'C' : return 1; // C_Tail
            case 'G' : return 2; // G_Tail
            case 'U' : return 3; // U_Tail
        }
    }

    static std::string seqT2U( std::string seq )
    {
        for( auto& nt : seq ) switch( nt )
            case 'T' : nt = 'U';
        return seq;
    }

    static std::string find_abundance( std::map< std::string, std::map< std::size_t, std::map< std::size_t, double >>>& seeds )
    {
        double sum = 0;
        double max_count = 0;
        std::string abundance = "";

        for( auto& seed : seeds )
        {
            sum = 0;
            for( auto& tail : seed.second ) for( auto& len : tail.second )
            {
                sum += len.second;
            }

            if( sum > max_count )
            {
                abundance = seed.first;
                max_count = sum;
            }
        }

        return abundance;
    }
    
    static std::set< std::pair< std::string, std::string >> make_anno_trimming(
            std::vector< ago::format::MDRawBed >& annotations,
            std::map< std::string, std::string >& genome_table,
            const std::size_t& filter_ppm,
            const std::string& biotype
            )
    {
        std::set< std::pair< std::string, std::string >> anno_trimmed;
        std::pair< std::string, std::string > anno_pair;

        std::map< std::pair< std::string, std::string >, double > anno_map;
        std::vector< std::pair< std::pair< std::string, std::string >, double >> anno_vec;

        std::string anno_first;
        std::string anno_second;
        bool isbiotype = false;

        std::size_t trim_top;
        std::size_t trim_last;

        for( auto& raw_bed : annotations )
        {
            if( raw_bed.ppm_ < filter_ppm ) continue;
            isbiotype = ( biotype == "" ? true : false );
            if( 0 < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[0].empty() )
            {
                if( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[0][0] == "miRNA" || raw_bed.annotation_info_[0][0] == "mirtron" ))
                    isbiotype = true;

                if( biotype != "" && raw_bed.annotation_info_[0][0] == biotype )
                    isbiotype = true;

                if( !isbiotype ) continue;

                for( std::size_t j = 0; j < raw_bed.annotation_info_[0].size(); j+=2 )
                {
                    anno_first  = raw_bed.annotation_info_[0][ j+1 ];
                    anno_second = GeneTypeAnalyzerCounting::seqT2U( raw_bed.getReadSeq( genome_table ).substr( 1, 7 ))
                                + ( GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag ) != "" ? ( "|" + GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag )) : "" );
                    anno_pair   = std::make_pair( anno_first, anno_second );

                    if( anno_map.find( anno_pair ) == anno_map.end() ) anno_map[ anno_pair ] = 0.0;
                    anno_map[ anno_pair ] += raw_bed.ppm_;
                }
            }
        }

        for( auto& anno : anno_map ) anno_vec.emplace_back( anno );

        std::sort( anno_vec.begin(), anno_vec.end(),
            []( const std::pair< std::pair< std::string, std::string >, double >& a,
                const std::pair< std::pair< std::string, std::string >, double >& b )
            {
                if( a.second == b.second )
                    return a.first > b.first;
                else
                    return a.second > b.second;
            });

        trim_top  = anno_vec.size() * 15 / 100;
        trim_last = anno_vec.size() * 85 / 100;

        for( std::size_t i = trim_top -1; i < trim_last; ++i )
            anno_trimmed.emplace( anno_vec[i].first );

        return anno_trimmed;
    }

    static void get_trimming( auto& anno_table_tail, auto& anno_table_trim, const double& percentage )
    {
        std::size_t trim_top;
        std::size_t trim_last;
        double smallest = 9999.99;

        std::set< std::string > trims_set;
        std::map< std::string, double > trims_map;
        std::vector< std::pair< std::string, double >> trims_vec;

        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            for( std::size_t i = 0; i < anno_table_tail[ smp ].size(); ++i )
            {
                for( auto& anno : anno_table_tail[ smp ][i] )
                {
                    if( trims_map.find( anno.first ) == trims_map.end() ) trims_map[ anno.first ] = 0.0;
                    for( auto& len : anno.second ) trims_map[ anno.first ] += len.second;
                }
            }

            for( auto& anno : trims_map )
                if( anno.second < smallest )
                    smallest = anno.second;

            for( auto& anno : trims_map )
                if( anno.second > smallest )
                    trims_vec.emplace_back( anno );

            smallest = 9999.99;
            trims_map.clear();

            std::sort( trims_vec.begin(), trims_vec.end(),
                []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
                {
                    if( a.second == b.second )
                        return a.first > b.first;
                    else
                        return a.second > b.second;
                });

            trim_top  = trims_vec.size() * 15 / 100;
            trim_last = trims_vec.size() * 85 / 100;

            for( std::size_t i = trim_top -1; i < trim_last; ++i )
                trims_set.emplace( trims_vec[i].first );

            trims_vec.clear();

            for( std::size_t i = 0; i < anno_table_tail[ smp ].size(); ++i )
            {
                for( auto& anno : anno_table_tail[ smp ][i] )
                {
                    if( trims_set.find( anno.first ) == trims_set.end() ) continue;

                    for( auto& len : anno.second )
                    {
                        if( anno_table_trim[ smp ][ i ][ anno.first ].find( len.first ) == anno_table_trim[ smp ][ i ][ anno.first ].end() )
                            anno_table_trim[ smp ][ i ][ anno.first ][ len.first ] = 0.0;
                            anno_table_trim[ smp ][ i ][ anno.first ][ len.first ] += len.second;
                    }
                }
            }

            trims_set.clear();
        }
    }

    static void make_anno_table(
            // std::vector< AnnotationRawBed<> >& annotations,
            std::vector< ago::format::MDRawBed >& annotations,
            std::vector< CountingTableType >& anno_table,
            std::map< std::string, std::string >& anno_mark,
            std::map< std::string, std::string >& genome_table,
            const std::size_t filter_ppm,
            const std::string biotype = "", // default annotation type
            const bool trimming = false
            )
    {
        //                  annotation      seed          count
        std::map< std::pair< std::string, std::string >, double > anno_check;
        std::map< std::string, std::map< std::string, std::map< std::size_t, std::map< std::size_t, double >>>> anno_temp;
        //          annotation              seed                    tail    count
        
        std::pair< std::string, std::string > anno_pair;

        std::string anno_idx;
        std::string anno_first;
        std::string anno_second;

        std::string gene_name;
        std::string gene_seed;

        std::string abundance;
        std::size_t read_len;
        std::size_t tail;

        double anno_counts;
        bool isbiotype = false;

        std::set< std::pair< std::string, std::string >> anno_trimmed;
        if( trimming ) anno_trimmed = make_anno_trimming( annotations, genome_table, filter_ppm, biotype );

        for( auto& raw_bed : annotations )
        {
            if( raw_bed.ppm_ < filter_ppm ) continue;
            isbiotype = ( biotype == "" ? true : false );
            read_len = raw_bed.length_ - raw_bed.tail_length_;
            tail = which_tail( GeneTypeAnalyzerCounting::seqT2U( raw_bed.getTail() ));

            // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
            {
                std::size_t i = 0; // do first priority
                if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                {
                    if( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))
                        isbiotype = true;

                    if( biotype != "" && raw_bed.annotation_info_[i][0] == biotype )
                        isbiotype = true;

                    if( !isbiotype ) continue;

                    for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                    {
                        gene_name = raw_bed.annotation_info_[i][ j+1 ];
                        gene_seed = GeneTypeAnalyzerCounting::seqT2U( raw_bed.getReadSeq( genome_table ).substr( 1, 7 ))
                                  + ( GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag ) != "" ? ( "|" + GeneTypeAnalyzerCounting::seqT2U( raw_bed.seed_md_tag )) : "" );

                        if( trimming && anno_trimmed.find( std::make_pair( gene_name, gene_seed )) == anno_trimmed.end() )
                            continue;

                        anno_first  = biotype == "" ? raw_bed.annotation_info_[i][j] : gene_name;
                        anno_second = biotype == "" ? "" : gene_seed;

                        anno_pair   = std::make_pair( anno_first, anno_second );

                        if( biotype != "" && raw_bed.is_filtered_ != 0 ) anno_mark[ anno_first + "_" + anno_second ] = "!";
                        
                        if( anno_temp[ anno_first ][ anno_second ][ tail ].find( read_len ) ==
                            anno_temp[ anno_first ][ anno_second ][ tail ].end() )
                            anno_temp[ anno_first ][ anno_second ][ tail ][ read_len ] = 0.0;

                        anno_counts = raw_bed.ppm_;
                        // anno_counts = biotype == "" ? raw_bed.ppm_ : ((double)(raw_bed.reads_count_)/(double)(raw_bed.multiple_alignment_site_count_));

                        if( anno_check.find( anno_pair ) == anno_check.end() )
                            anno_check[ anno_pair ] = 0.0;

                        anno_check[ anno_pair ] += anno_counts;
                    }
                }
            }

            for( auto& anno : anno_check )
            {
                anno_temp[ anno.first.first ][ anno.first.second ][ tail ][ read_len ] += anno.second;
            }

            anno_check.clear();
        }

        for( auto& anno : anno_temp )
        {
            abundance = find_abundance( anno.second );
            if( biotype != "" ) anno_mark[ anno.first + "_" + abundance ] += "*";

            for( auto& seed : anno.second )
            {
                anno_idx = biotype == "" ? anno.first : ( anno.first + "_" + seed.first );

                for( auto& tal : seed.second ) for( auto& len : tal.second )
                {
                    if( anno_table[ tal.first ][ anno_idx ].find( len.first ) == anno_table[ tal.first ][ anno_idx ].end() )
                        anno_table[ tal.first ][ anno_idx ][ len.first ] = 0.0;

                    anno_table[ tal.first ][ anno_idx ][ len.first ] += len.second;
                }
            }
        }
    }

    static void make_ana_table(
            std::vector< CountingTableType >& anno_table,
            std::map< std::string, std::string >& anno_mark,
            std::map< std::string, std::map< std::string, double >>& seed_match_table,
            const std::string& token
            )
    {
        std::vector< std::string > split;
        std::vector< CountingTableType > anno_table_temp = std::vector< CountingTableType >( anno_table.size() );
        std::map< std::string, std::string > anno_mark_temp;

        bool is_seed = token == "Seed" ? true : false;
        bool is_arms = token == "3p" || token == "5p" ? true : false;
        bool is_lens = !is_seed && !is_arms ? true : false;

        for( std::size_t tail = 0; tail < anno_table.size(); tail++ )
        {
            for( auto& name_seed : anno_table[ tail ] )
            {
                boost::iter_split( split, name_seed.first, boost::algorithm::first_finder( "_" ));

                if( split.size() > 2 )
                {
                    for( std::size_t i = 1; i < split.size() -2; ++i )
                        split[0] = split[0] + "_" + split[i];
                    split[1] = split[ split.size() -1 ];
                }

                if( !is_seed )
                {
                    if( is_arms && get_arm( split[0] ) != token )
                        continue;

                    for( auto& len : name_seed.second )
                    {
                        if( is_lens && token != std::to_string( len.first ))
                            continue;

                        if( anno_table_temp[ tail ][ name_seed.first ].find( len.first ) == anno_table_temp[ tail ][ name_seed.first ].end() )
                            anno_table_temp[ tail ][ name_seed.first ][ len.first ] = 0.0;

                        anno_table_temp[ tail ][ name_seed.first ][ len.first ] += len.second;
                    }

                    if( anno_mark.find( name_seed.first ) != anno_mark.end() )
                        anno_mark_temp[ name_seed.first ] = anno_mark[ name_seed.first ];

                    continue;
                }

                for( auto& len : name_seed.second )
                {
                    if( anno_table_temp[ tail ][ split[1] ].find( len.first ) == anno_table_temp[ tail ][ split[1] ].end() )
                        anno_table_temp[ tail ][ split[1] ][ len.first ] = 0.0;

                    anno_table_temp[ tail ][ split[1] ][ len.first ] += len.second;

                    if( seed_match_table[ split[1] ].find( split[0] ) == seed_match_table[ split[1] ].end() )
                        seed_match_table[ split[1] ][ split[0] ] = 0.0;

                    seed_match_table[ split[1] ][ split[0] ] += len.second;
                }
            }
        }

        if( !is_seed )
            anno_mark = std::move( anno_mark_temp );

        anno_table = std::move( anno_table_temp );
    }

    static std::string get_arm( const std::string& anno )
    {
        std::vector< std::string > split;
        boost::iter_split( split, anno, boost::algorithm::first_finder( "-" ));
        return split[ split.size() -1 ];
    }

    static void table_refinding(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::size_t& min_len,
            const std::size_t& max_len,
            const double& sudo_count = 0.0
            )
    {
        std::set< std::size_t > len_idx;

        std::vector< std::vector< CountingTableType >> anno_table_tail_temp
            ( anno_table_tail.size(), std::vector< CountingTableType >( anno_table_tail[0].size() ));

        for( auto& len : ano_len_idx.second )
        {
            if( min_len != 0 && len < min_len ) continue;
            if( max_len != 0 && len > max_len ) continue;

            len_idx.emplace( len );
        }

        if( min_len != 0 || max_len != 0 ) ano_len_idx.second = len_idx;

        for( auto& anno : ano_len_idx.first )
        {
            for( auto& len : ano_len_idx.second )
            {
                if( min_len != 0 && len < min_len ) continue;
                if( max_len != 0 && len > max_len ) continue;

                for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
                {
                    for( std::size_t tal = 0; tal < anno_table_tail[ smp ].size(); ++tal )
                    {
                        if( min_len == 0 && max_len == 0 )
                        {
                            if( anno_table_tail[ smp ][ tal ].find( anno )        == anno_table_tail[ smp ][ tal ].end() ||
                                anno_table_tail[ smp ][ tal ][ anno ].find( len ) == anno_table_tail[ smp ][ tal ][ anno ].end() )
                                anno_table_tail[ smp ][ tal ][ anno ][ len ] = sudo_count;

                            continue;
                        }

                        if( anno_table_tail_temp[ smp ][ tal ][ anno ].find( len ) == anno_table_tail_temp[ smp ][ tal ][ anno ].end() )
                            anno_table_tail_temp[ smp ][ tal ][ anno ][ len ] = sudo_count;

                        if( anno_table_tail[ smp ][ tal ][ anno ].find( len ) != anno_table_tail[ smp ][ tal ][ anno ].end() )
                            anno_table_tail_temp[ smp ][ tal ][ anno ][ len ]  = anno_table_tail[ smp ][ tal ][ anno ][ len ];
                    }
                }
            }
        }

        if( min_len != 0 || max_len != 0 ) anno_table_tail = anno_table_tail_temp;
    }

    static void seed_refinding(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::map< std::string, double >>>& seed_match_table
            )
    {
        std::set< std::string > anno_list;
        std::vector< std::map< std::string, std::map< std::string, double >>> seed_match_temp
            ( anno_table_tail.size(), std::map< std::string, std::map< std::string, double >>() );

        for( auto& seed : ano_len_idx.first )
            for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
                if( seed_match_table[ smp ].find( seed ) != seed_match_table[ smp ].end() )
                    for( auto& anno : seed_match_table[ smp ][ seed ] )
                        anno_list.emplace( anno.first );

        for( auto& seed : ano_len_idx.first )
        {
            for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp ) for( auto& anno : anno_list )
            {
                seed_match_temp[ smp ][ seed ][ anno ] = 
                    seed_match_table[ smp ][ seed ].find( anno ) == seed_match_table[ smp ][ seed ].end()
                    ? 0.0 : seed_match_table[ smp ][ seed ][ anno ];
            }
        }

        seed_match_table = std::move( seed_match_temp );
    }

    static auto get_quantiled_ppm_map( std::vector< std::vector< CountingTableType >>& anno_table_tail )
    {
        std::vector< std::pair<              // each samples
            std::map< std::string, double >, // mir
            std::map< std::string, double >  // iosmir
                >> quantiled_ppm( anno_table_tail.size() );

        std::vector< std::string > split;

        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            for( std::size_t i = 0; i < anno_table_tail[ smp ].size(); ++i )
            {
                for( auto& anno : anno_table_tail[ smp ][i] )
                {
                    boost::iter_split( split, anno.first, boost::algorithm::first_finder( "_" ));

                    for( std::size_t i = 1; i < split.size() -2; ++i )
                        split[0] = split[0] + "_" + split[i];

                    if( quantiled_ppm[ smp ].first.find( split[0] ) == quantiled_ppm[ smp ].first.end() )
                        quantiled_ppm[ smp ].first[ split[0] ] = 0.0;

                    if( quantiled_ppm[ smp ].second.find( anno.first ) == quantiled_ppm[ smp ].second.end() )
                        quantiled_ppm[ smp ].second[ anno.first ] = 0.0;

                    for( auto& len : anno.second )
                    {
                        quantiled_ppm[ smp ].first[ split[0] ] += len.second;
                        quantiled_ppm[ smp ].second[ anno.first ] += len.second;
                    }
                }
            }
        }

        return quantiled_ppm;
    }
};

} // end of namespace algorithm
} // end of namespace ago
