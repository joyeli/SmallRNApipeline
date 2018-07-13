#pragma once

namespace ago {
namespace algorithm {

struct SA_Type
{
    std::vector< double > GMPM; // GMPM, GM, PM
    std::vector< double > Tail; // Atail, Ctail, Gtail, Other
    std::vector< double > tail; // aTail, cTail, gTail, Other is splite in to each Tail
    std::vector< std::map< char, double >> Seed; // 5'End-NT + Seed + 3'End-NT ( 1 + 7 + 1 )

    SA_Type()
        : GMPM( 3, 0.0 )
        , Tail( 5, 0.0 )
        , tail( 4, 0.0 )
        , Seed( 9, std::map< char, double >() )
    {}
};

class GeneTypeAnalyzerSA_Plot
{
  public:

    static std::vector< std::map< std::string, SA_Type >> anno_sa_table;

    GeneTypeAnalyzerSA_Plot()
    {}

    static void make_sa_plot_table(
            const std::string& biotype,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< BedSampleType >& bed_samples,
            auto& genome_table,
            const bool& isSeed = false
            )
    {
        anno_sa_table.clear();
        std::map< std::string, SA_Type > anno_map;

        for( auto& anno : ano_len_idx.first )
            anno_map.emplace( anno, SA_Type() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            anno_sa_table.emplace_back( anno_map );
            make_sa_counting_table( biotype, bed_samples[ smp ], anno_sa_table[ smp ], genome_table, isSeed );
            // debug( anno_sa_table[ smp ] );
            formation_sa_table( anno_sa_table[ smp ] );
        }
    }

    static bool check_biotype( const auto& raw_bed, const auto& biotype )
    {
        if( raw_bed.annotation_info_[0][0] == biotype || ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[0][0] == "miRNA" || raw_bed.annotation_info_[0][0] == "mirtron" )))
            return true;
        return false;
    }

    static void insert_sa_table(
            const auto& sa_it,
            const auto& sequence,
            const auto& tail_seq,
            const auto& tail,
            const auto& ppm,
            auto& md_map,
            auto& tc_set
            )
    {
        sa_it->second.GMPM[0] += ppm;

        if( tail != 5 )
        {
            sa_it->second.GMPM[2] += ppm;
            sa_it->second.Tail[ tail ] += ppm;

            if( tail == 4 )
            {
                std::map< char, std::size_t > tail_count;

                for( std::size_t j = 0; j < tail_seq.length(); ++j )
                {
                    if( tail_count.find( tail_seq.at(j) ) == tail_count.end() )
                        tail_count[ tail_seq.at(j) ] = 0;
                    tail_count[ tail_seq.at(j) ] += 1;
                }

                for( auto& tail_char : tail_count ) switch ( tail_char.first )
                {
                    case 'A' : sa_it->second.tail[0] += (double)(tail_char.second) / (double)(tail_seq.length()) * ppm; break;
                    case 'C' : sa_it->second.tail[1] += (double)(tail_char.second) / (double)(tail_seq.length()) * ppm; break;
                    case 'G' : sa_it->second.tail[2] += (double)(tail_char.second) / (double)(tail_seq.length()) * ppm; break;
                    case 'T' : sa_it->second.tail[3] += (double)(tail_char.second) / (double)(tail_seq.length()) * ppm; break;
                }
            }
            else sa_it->second.tail[ tail ] += ppm;
        }
        else sa_it->second.GMPM[1] += ppm;

        char nt = ' ';

        for( std::size_t j = 0; j < 8; ++j )
        {
            nt = sequence.at(j);

            if( md_map.find(j) != md_map.end() ) nt = md_map[j];
            if( tc_set.find(j) != tc_set.end() ) nt = 'T';

            if( sa_it->second.Seed[j].find( nt ) == sa_it->second.Seed[j].end() )
                sa_it->second.Seed[j][ nt ] = 0;

            sa_it->second.Seed[j][ nt ] += ppm;
        }

        nt = sequence.at( sequence.length() -1 );

        if( md_map.find( sequence.length() -1 ) != md_map.end() ) nt = md_map[ sequence.length() -1 ];
        if( tc_set.find( sequence.length() -1 ) != tc_set.end() ) nt = 'T';

        if( sa_it->second.Seed[8].find( nt ) == sa_it->second.Seed[8].end() )
            sa_it->second.Seed[8][ nt ] = 0;

        sa_it->second.Seed[8][ nt ] += ppm;
    }

    static void make_sa_counting_table(
            const std::string& biotype,
            BedSampleType& bed_sample,
            std::map< std::string, SA_Type >& sa_table,
            auto& genome_table,
            const bool& isSeed = false
            )
    {
        std::string gene_name;
        std::string gene_seed;

        std::string sequence;
        std::string tail_seq;

        std::size_t tail;

        for( auto& raw_bed : bed_sample.second )
        {
            if( !raw_bed.annotation_info_.empty() && !raw_bed.annotation_info_[0].empty() )
            {
                if( !check_biotype( raw_bed, biotype )) continue;

                sequence = raw_bed.getReadSeq( genome_table );
                tail_seq = raw_bed.getTail();

                tail = GeneTypeAnalyzerCounting::which_tail( tail_seq );

                for( std::size_t i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    gene_seed = raw_bed.getReadSeq( genome_table ).substr( 1, 7 )
                            + ( raw_bed.seed_md_tag != "" ? ( "|" + raw_bed.seed_md_tag ) : "" );

                    gene_name = isSeed ? gene_seed : ( raw_bed.annotation_info_[0][ i+1 ] + "_" + gene_seed );

                    if( sa_table.find( gene_name ) != sa_table.end() )
                        insert_sa_table( sa_table.find( gene_name ), sequence, tail_seq, tail, raw_bed.ppm_, raw_bed.md_map, raw_bed.tc_set );
                }
            }
        }
    }

    static SA_Type get_sum( std::map< std::string, SA_Type >& sa_table )
    {
        SA_Type sum;

        for( std::size_t i = 0; i < 9; ++i ) sum.Seed[i]['S'] = 0.0;

        for( auto& sa : sa_table )
        {
            for( std::size_t i = 0; i < 5; ++i ) sum.Tail[i] += sa.second.Tail[i];
            for( std::size_t i = 0; i < 4; ++i ) sum.tail[i] += sa.second.tail[i];
            for( std::size_t i = 0; i < 9; ++i )
            for( auto & nt : sa.second.Seed[i] ) sum.Seed[i]['S'] += nt.second;
        }

        return sum;
    }

    static void formation_sa_table( std::map< std::string, SA_Type >& sa_table )
    {
        // SA_Type sum = get_sum( sa_table );

        for( auto& sa : sa_table )
        {
            // for( std::size_t i = 0; i < 5; ++i ) sa.second.Tail[i] = sa.second.Tail[i] / sum.Tail[i];
            // for( std::size_t i = 0; i < 4; ++i ) sa.second.tail[i] = sa.second.tail[i] / sum.tail[i];
            for( std::size_t i = 0; i < 9; ++i )
            {
                // for( auto & nt : sa.second.Seed[i] ) nt.second = nt.second / sum.Seed[i]['S'];
                if( sa.second.Seed[i].find( 'A' ) == sa.second.Seed[i].end() ) sa.second.Seed[i][ 'A' ] = 0.0;
                if( sa.second.Seed[i].find( 'C' ) == sa.second.Seed[i].end() ) sa.second.Seed[i][ 'C' ] = 0.0;
                if( sa.second.Seed[i].find( 'G' ) == sa.second.Seed[i].end() ) sa.second.Seed[i][ 'G' ] = 0.0;
                if( sa.second.Seed[i].find( 'T' ) == sa.second.Seed[i].end() ) sa.second.Seed[i][ 'T' ] = 0.0;
            }
        }
    }

    static void output_sa_plot(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            const std::string& sample_name,
            const std::size_t& smp
            )
    {
        std::ofstream outtll( output_name + sample_name + "-tail.tsv" );
        std::ofstream outtlo( output_name + sample_name + "-tail_with_other.tsv" );
        std::ofstream outnt0( output_name + sample_name + "-nt_0.tsv" );
        std::ofstream outnt1( output_name + sample_name + "-nt_1.tsv" );
        std::ofstream outnt2( output_name + sample_name + "-nt_2.tsv" );
        std::ofstream outnt3( output_name + sample_name + "-nt_3.tsv" );
        std::ofstream outnt4( output_name + sample_name + "-nt_4.tsv" );
        std::ofstream outnt5( output_name + sample_name + "-nt_5.tsv" );
        std::ofstream outnt6( output_name + sample_name + "-nt_6.tsv" );
        std::ofstream outnt7( output_name + sample_name + "-nt_7.tsv" );
        std::ofstream outnt8( output_name + sample_name + "-nt_last.tsv" );

        outtll << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outtlo << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT\tO";
        outnt0 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt1 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt2 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt3 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt4 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt5 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt6 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt7 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
        outnt8 << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";

        for( auto& anno : ano_len_idx.first )
        {
            if( anno_sa_table[ smp ].find( anno ) == anno_sa_table[ smp ].end() ) continue;

            outtll << "\n" << anno;
            outtlo << "\n" << anno;
            outnt0 << "\n" << anno;
            outnt1 << "\n" << anno;
            outnt2 << "\n" << anno;
            outnt3 << "\n" << anno;
            outnt4 << "\n" << anno;
            outnt5 << "\n" << anno;
            outnt6 << "\n" << anno;
            outnt7 << "\n" << anno;
            outnt8 << "\n" << anno;

            for( std::size_t i = 0; i < 3; ++i )
            {
                outtll << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outtlo << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt0 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt1 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt2 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt3 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt4 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt5 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt6 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt7 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
                outnt8 << "\t" << anno_sa_table[ smp ][ anno ].GMPM[i];
            }

            for( std::size_t i = 0; i < 4; ++i ) outtll << "\t" << anno_sa_table[ smp ][ anno ].tail[i];
            for( std::size_t i = 0; i < 5; ++i ) outtlo << "\t" << anno_sa_table[ smp ][ anno ].Tail[i];
            for( std::size_t i = 0; i < 9; ++i ) switch( i )
            {
                case 0 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt0 << "\t" << nt.second; break;
                case 1 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt1 << "\t" << nt.second; break;
                case 2 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt2 << "\t" << nt.second; break;
                case 3 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt3 << "\t" << nt.second; break;
                case 4 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt4 << "\t" << nt.second; break;
                case 5 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt5 << "\t" << nt.second; break;
                case 6 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt6 << "\t" << nt.second; break;
                case 7 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt7 << "\t" << nt.second; break;
                case 8 : for( auto& nt : anno_sa_table[ smp ][ anno ].Seed[i] ) outnt8 << "\t" << nt.second; break;
            }
        }

        outtll.close();
        outtlo.close();
        outnt0.close();
        outnt1.close();
        outnt2.close();
        outnt3.close();
        outnt4.close();
        outnt5.close();
        outnt6.close();
        outnt7.close();
        outnt8.close();
    }

    static void debug( std::map< std::string, SA_Type >& sa_table )
    {
        for( auto& anno : sa_table )
        {
            auto& sa = anno.second;
            std::cerr << anno.first;

            for( std::size_t i = 0; i < 3; ++i ) std::cerr << "\t" << sa.GMPM[i];
            for( std::size_t i = 0; i < 5; ++i ) std::cerr << "\t" << sa.Tail[i];
            for( std::size_t i = 0; i < 4; ++i ) std::cerr << "\t" << sa.tail[i];
            for( std::size_t i = 0; i < 9; ++i )
            {
                std::cerr << ( sa.Seed[i].find( 'A' ) == sa.Seed[i].end() ? "\tAnotFound" : ( "\t" + std::to_string( sa.Seed[i][ 'A' ] )));
                std::cerr << ( sa.Seed[i].find( 'C' ) == sa.Seed[i].end() ? "\tCnotFound" : ( "\t" + std::to_string( sa.Seed[i][ 'C' ] )));
                std::cerr << ( sa.Seed[i].find( 'G' ) == sa.Seed[i].end() ? "\tGnotFound" : ( "\t" + std::to_string( sa.Seed[i][ 'G' ] )));
                std::cerr << ( sa.Seed[i].find( 'T' ) == sa.Seed[i].end() ? "\tTnotFound" : ( "\t" + std::to_string( sa.Seed[i][ 'T' ] )));
            }

            std::cerr << "\n";
        }
    }

    static void output_sa_plot_visualization( const std::string& output_name, const std::string& biotype, const bool& isSeed )
    {
    }

};

std::vector< std::map< std::string, SA_Type >> GeneTypeAnalyzerSA_Plot::anno_sa_table;

} // end of namespace algorithm
} // end of namespace ago