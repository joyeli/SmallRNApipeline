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
            const std::string& biotype = "" // default annotation type
            )
    {
        std::vector< std::set< std::string >> smp_ano_idx( bed_samples.size(), std::set< std::string >() );
        std::vector< std::set< std::size_t >> smp_len_idx( bed_samples.size(), std::set< std::size_t >() );

        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &genome_table, &biotype, &bed_samples, &smp_ano_idx, &smp_len_idx ] () mutable
            {
                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    smp_len_idx[ smp ].emplace( raw_bed.length_ - raw_bed.tail_length_ );

                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                        {
                            if( biotype != "" && raw_bed_info[i] != biotype ) continue;
                            smp_ano_idx[ smp ].emplace( raw_bed_info[ i + ( biotype == "" ? 0 : 1 )]
                                + ( biotype == "" ? "" :
                                (
                                    "_" + raw_bed.getReadSeq( genome_table ).substr( 1, 7 )
                                    + ( raw_bed.seed_md_tag != "" ? ( "|" + raw_bed.seed_md_tag ) : "" )
                                ) ));
                        }
                    }
                }
            });
        }

        smp_parallel_pool.flush_pool();

        std::set< std::string > ano_idx;
        std::set< std::size_t > len_idx;
        std::set< std::size_t > len_temp;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            ano_idx.insert( smp_ano_idx[ smp ].begin(), smp_ano_idx[ smp ].end() );
            len_idx.insert( smp_len_idx[ smp ].begin(), smp_len_idx[ smp ].end() );
        }

        std::size_t count = 0;
        std::size_t last_len;

        for( auto& len : len_idx )
        {
            if( count == 0 ) last_len = len -1;
            if( len - last_len > 1 )
                for( std::size_t i = last_len + 1; i < len; ++i )
                    len_temp.emplace( i ); 
            count++;
        }

        len_idx.insert( len_temp.begin(), len_temp.end() );

        return std::make_pair( ano_idx, len_idx );
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
            case 'T' : return 3; // T_Tail
        }
    }

    // static double get_ppm( std::vector< AnnotationRawBed<> >& annotations, double ppm = 1000000 )
    static double get_ppm( std::vector< ago::format::MDRawBed >& annotations, double ppm = 1000000 )
    {
        double sum = 0;

        for( auto& anno : annotations )
        {
            for( auto& info : anno.annotation_info_ )
            {
                for( int i = 0; i < info.size(); i+=2 )
                {
                    sum += ( anno.reads_count_ / anno.multiple_alignment_site_count_ );
                }
            }
        }

        return ppm / sum;
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

    static void make_anno_table(
            // std::vector< AnnotationRawBed<> >& annotations,
            std::vector< ago::format::MDRawBed >& annotations,
            std::vector< CountingTableType >& anno_table,
            std::map< std::string, std::string >& anno_mark,
            std::map< std::string, std::string >& genome_table,
            const std::string& biotype = "" // default annotation type
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
        std::string abundance;

        std::size_t tail;
        std::size_t read_len;

        double anno_check_sum;
        double anno_counts;
        double ppm = get_ppm( annotations );

        for( auto& raw_bed : annotations )
        {
            anno_check.clear();
            anno_check_sum = 0.0;

            read_len = raw_bed.length_ - raw_bed.tail_length_;
            tail = which_tail( raw_bed.getTail() );

            for( auto& raw_bed_info : raw_bed.annotation_info_ )
            {
                for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                {
                    if( biotype != "" && raw_bed_info[i] != biotype ) continue;

                    anno_first  = biotype == "" ? raw_bed_info[i] : raw_bed_info[ i+1 ];
                    anno_second = biotype == "" ? ""
                                :   raw_bed.getReadSeq( genome_table ).substr( 1, 7 )
                                + ( raw_bed.seed_md_tag != "" ? ( "|" + raw_bed.seed_md_tag ) : "" );
                    anno_pair   = std::make_pair( anno_first, anno_second );

                    if( biotype != "" && raw_bed.is_filtered_ != 0 ) anno_mark[ anno_first + "_" + anno_second ] = "!";
                    
                    if( anno_temp[ anno_first ][ anno_second ][ tail ].find( read_len ) ==
                        anno_temp[ anno_first ][ anno_second ][ tail ].end() )
                        anno_temp[ anno_first ][ anno_second ][ tail ][ read_len ] = 0.0;

                    anno_check_sum++;
                    anno_counts = raw_bed.reads_count_ * ppm / raw_bed.multiple_alignment_site_count_;

                    if( anno_check.find( anno_pair ) == anno_check.end() ) anno_check[ anno_pair ] = 0;
                    anno_check[ anno_pair ] += anno_counts;
                }
            }

            for( auto& anno : anno_check )
            {
                anno.second = anno.second / anno_check_sum;
                anno_temp[ anno.first.first ][ anno.first.second ][ tail ][ read_len ] += anno.second;
            }
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
};

} // end of namespace algorithm
} // end of namespace ago
