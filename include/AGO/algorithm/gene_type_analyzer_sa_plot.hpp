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

    static void output_sa_plot_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $Type = $_POST['Type'];" << "\n";
        output << "        $isBin = $_POST['isBin'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $Ranking = $_POST['Ranking'];" << "\n";
        output << "        $nNumber = $_POST['nNumber'];" << "\n";
        output << "        $isRatio = $_POST['isRatio'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $isTrimmed = $_POST['isTrimmed'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Type ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Type onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Type=='') echo 'selected'; echo 'value= >Select Type</option>';" << "\n";
        output << "" << "\n";
        output << "        $Type_List = array( '5’End', '3’End', 'Seed0', 'Seed1', 'Seed2', 'Seed3', 'Seed4', 'Seed5', 'Seed6', 'Tail', 'TailwithOther' );" << "\n";
        output << "        $Type_Size = Count( $Type_List );" << "\n";
        output << "        $Type_Name = '';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Type_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Type_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Type == $Type_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Type_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        switch( $Type )" << "\n";
        output << "        {" << "\n";
        output << "            case '5’End'         : $Type_Name = '-nt_0.tsv'           ; break;" << "\n";
        output << "            case '3’End'         : $Type_Name = '-nt_last.tsv'        ; break;" << "\n";
        output << "            case 'Seed0'         : $Type_Name = '-nt_1.tsv'           ; break;" << "\n";
        output << "            case 'Seed1'         : $Type_Name = '-nt_2.tsv'           ; break;" << "\n";
        output << "            case 'Seed2'         : $Type_Name = '-nt_3.tsv'           ; break;" << "\n";
        output << "            case 'Seed3'         : $Type_Name = '-nt_4.tsv'           ; break;" << "\n";
        output << "            case 'Seed4'         : $Type_Name = '-nt_5.tsv'           ; break;" << "\n";
        output << "            case 'Seed5'         : $Type_Name = '-nt_6.tsv'           ; break;" << "\n";
        output << "            case 'Seed6'         : $Type_Name = '-nt_7.tsv'           ; break;" << "\n";
        output << "            case 'Tail'          : $Type_Name = '-tail.tsv'           ; break;" << "\n";
        output << "            case 'TailwithOther' : $Type_Name = '-tail_with_other.tsv'; break;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $Sample_List = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_List ) -1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Sample = Explode( '-', $TSV_List[$i] );" << "\n";
        output << "            $Sample_List[ $Sample[0] ] = $Sample[0];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select Tsv</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Sample_List as $tsv => $smp )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$tsv.' ';" << "\n";
        output << "            if( $TSV_File == $tsv ) echo 'selected ';" << "\n";
        output << "            echo '>'.$tsv.'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "" << "\n";
        output << "        $nCount = 0;" << "\n";
        output << "        $Header = Array();" << "\n";
        output << "" << "\n";
        output << "        $Index_Array = Array();" << "\n";
        output << "        $Value_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        if( File_Exists( $TSV_File.$Type_Name ))" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = new SplFileObject( $TSV_File.$Type_Name );" << "\n";
        output << "            $isHeader = true;" << "\n";
        output << "" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Ranking == 'Ranking' || $isRatio == 'isRatio' ||" << "\n";
        output << "                    $Ranking == ''        || $isRatio == '' ) break;" << "\n";
        output << "" << "\n";
        output << "                $Ratio_Array = Array();" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                $Ranking_Type = 0;" << "\n";
        output << "" << "\n";
        output << "                if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "                $inFile_Lines = Rtrim( $inFile_Lines );" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "                if( $isHeader )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Ranking == 'Tailing％' )" << "\n";
        output << "                        $Header[0] = $inFile_Line[2];" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 4; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                        $Header[ $i - ( $Ranking != 'Tailing％' ? 4 : 3 )] = $inFile_Line[$i];" << "\n";
        output << "" << "\n";
        output << "                    $isHeader = false;" << "\n";
        output << "                    continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Filter != '' && $Filter != 'FilterGMPM' && $inFile_Line[1] < $Filter  ) continue;" << "\n";
        output << "                if( $Ranking == 'Tailing％' && $inFile_Line[3] == 0 ) continue;" << "\n";
        output << "                if( $inFile_Line[1] == 0 ) continue;" << "\n";
        output << "" << "\n";
        output << "                switch( $Ranking )" << "\n";
        output << "                {" << "\n";
        output << "                    case 'GM'        : $Ranking_Type = $inFile_Line[2] ; break;" << "\n";
        output << "                    case 'PM'        : $Ranking_Type = $inFile_Line[3] ; break;" << "\n";
        output << "                    case 'Tailing％' : $Ranking_Type = $inFile_Line[3] / $inFile_Line[1]; break;" << "\n";
        output << "                    default          : $Ranking_Type = $inFile_Line[1] ; break;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Ranking == 'Tailing％' )" << "\n";
        output << "                    $Ratio_Array[0] = $inFile_Line[2];" << "\n";
        output << "" << "\n";
        output << "                For( $i = 4; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    $Ratio_Array[ $i - ( $Ranking != 'Tailing％' ? 4 : 3 )] = $inFile_Line[$i];" << "\n";
        output << "" << "\n";
        output << "                Array_Push( $Index_Array, $Ranking_Type );" << "\n";
        output << "                Array_Push( $Value_Array, $Ratio_Array );" << "\n";
        output << "                $nCount++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Array_Multisort( $Index_Array, $Value_Array );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Alln = Count( $Value_Array );" << "\n";
        output << "        $Trimming = $Alln * ( $isTrimmed == '' || $isTrimmed == 'isTrimmed' ? 0 : $isTrimmed ) / 100;" << "\n";
        output << "" << "\n";
        output << "#<!--================== Ranking ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Ranking onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Ranking=='') echo 'selected'; echo 'value= >Select Ranking</option>';" << "\n";
        output << "" << "\n";
        output << "        $Rank_List = array( 'GMPM', 'GM', 'PM', 'Tailing％' );" << "\n";
        output << "        $Rank_Size = Count( $Rank_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Rank_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Rank_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Ranking == $Rank_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Rank_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================= isTrimmed ===================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isTrimmed onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isTrimmed=='') echo 'selected'; echo 'value= >isTrimmed</option>';" << "\n";
        output << "" << "\n";
        output << "        $Trim_List = array( '5', '10', '15', '20', '25' );" << "\n";
        output << "        $Trim_Size = Count( $Trim_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Trim_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Trim_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isTrimmed == $Trim_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Trim_List[$i].'％</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Bin =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isBin onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isBin=='') echo 'selected'; echo 'value= >isBin</option>';" << "\n";
        output << "" << "\n";
        output << "        $isBin_List = array( 5, 10, 20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000 );" << "\n";
        output << "        $isBin_Size = Count( $isBin_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $isBin_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$isBin_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isBin == $isBin_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$isBin_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================ isRatio ==================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isRatio onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isRatio=='') echo 'selected'; echo 'value= >isRatio</option>';" << "\n";
        output << "" << "\n";
        output << "        $isRatio_List = array( 'Yes', 'No' );" << "\n";
        output << "        $isRatio_Size = Count( $isRatio_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $isRatio_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$isRatio_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isRatio == $isRatio_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$isRatio_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== nNumber =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        if( $nNumber == '' && $nCount == 0 )" << "\n";
        output << "            $nNumber = 'SetN';" << "\n";
        output << "" << "\n";
        output << "        if(( $nNumber == '' || $nNumber == 'Set N' || $nNumber > $nCount ) && $nCount != 0 )" << "\n";
        output << "            $nNumber = $nCount/2;" << "\n";
        output << "" << "\n";
        output << "        if( $isBin != '' && $isBin != 'isBin' )" << "\n";
        output << "            $nNumber = Round( $nNumber / $isBin ) * $isBin;" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=nNumber size=3 value='.$nNumber.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Filter =====================-->" << "\n";
        output << "        " << "\n";
        output << "        if( $Filter == '' ) $Filter = 'FilterGMPM';" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value='.$Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MakeTemp ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Sum_Array = Array();" << "\n";
        output << "        $Temp_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $L_Array = Array();" << "\n";
        output << "        $M_Array = Array();" << "\n";
        output << "        $R_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Header ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Sum_Array[$i] = 0.0;" << "\n";
        output << "            $L_Array[ $Header[$i] ] = Array();" << "\n";
        output << "            $M_Array[ $Header[$i] ] = Array();" << "\n";
        output << "            $R_Array[ $Header[$i] ] = Array();" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $AllArray = $Sum_Array;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $nNumber; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            Array_Push( $Temp_Array, $Value_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "            if( $isBin == '' || $isBin == 'isBin' || ( $i % $isBin == 0 && $i != 0 ) || $i == $nNumber -1 )" << "\n";
        output << "            {" << "\n";
        output << "                $N = ( $i == $nNumber -1 ? ( $i + 1 ) : $i );" << "\n";
        output << "                $Sum = 0.0;" << "\n";
        output << "                $SumArray = $Sum_Array;" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $Temp_Array ); ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $PM = 0.0;" << "\n";
        output << "" << "\n";
        output << "                    if( $isRatio == 'Yes' ) For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $L_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
        output << "" << "\n";
        output << "                $Temp_Array = Array();" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Alln; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $isSkip = false;" << "\n";
        output << "" << "\n";
        output << "            if(( $isTrimmed != '' && $isTrimmed != 'isTrimmed' ) &&" << "\n";
        output << "               ( $i < $Trimming || $i > $Alln - $Trimming )) $isSkip = true;" << "\n";
        output << "" << "\n";
        output << "            if( !$isSkip ) Array_Push( $Temp_Array, $Value_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "            if( $i == $Alln -1 )" << "\n";
        output << "            {" << "\n";
        output << "                $N = Count( $Temp_Array );" << "\n";
        output << "                $Sum = 0.0;" << "\n";
        output << "                $SumArray = $Sum_Array;" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $Temp_Array ); ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $PM = 0.0;" << "\n";
        output << "" << "\n";
        output << "                    if( $isRatio == 'Yes' ) For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $M_Array[ $Header[$j] ][0]  = $SumArray[$j];" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $M_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
        output << "" << "\n";
        output << "                $Temp_Array = Array();" << "\n";
        output << "                $Trimmedn = $N;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Label_Witdh = $Trimming == 0 ? 2.3 : 5.3;" << "\n";
        output << "        $Label = 'All = '.$Alln.( $Trimming == 0 ? '' : ( ', Trimmed = '.$Trimmedn ));" << "\n";
        output << "" << "\n";
        output << "        For( $n = $nNumber, $i = $Alln - $nNumber; $i < $Alln; $n--,++$i )" << "\n";
        output << "        {" << "\n";
        output << "            Array_Push( $Temp_Array, $Value_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "            if( $isBin == '' || $isBin == 'isBin' || ( $n % $isBin == 0 && $n != $nNumber ) || $i == $Alln -1 )" << "\n";
        output << "            {" << "\n";
        output << "                $N = ( $i == $Alln -1 ? ( $n - 1 ) : $n ) + $isBin;" << "\n";
        output << "                $Sum = 0.0;" << "\n";
        output << "                $SumArray = $Sum_Array;" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $Temp_Array ); ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $PM = 0.0;" << "\n";
        output << "" << "\n";
        output << "                    if( $isRatio == 'Yes' ) For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = 0; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = 0; $j < Count( $SumArray ); ++$j ) $R_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
        output << "" << "\n";
        output << "                $Temp_Array = Array();" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== FormateTemp ==================-->" << "\n";
        output << "" << "\n";
        output << "        Foreach( $L_Array as $ACGT => $N_Array )" << "\n";
        output << "        {" << "\n";
        output << "            $i = 0;" << "\n";
        output << "            $n = Floor( $nNumber / ( Count( $N_Array ) -1 ));" << "\n";
        output << "            $l = $n * ( Count( $N_Array ) -1 );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $N_Array as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                Unset( $L_Array[ $ACGT ][$N] );" << "\n";
        output << "                $L_Array[ $ACGT ][ $i != $l ? $i : $nNumber ] = $Value;" << "\n";
        output << "                $i += $n;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Foreach( $R_Array as $ACGT => $N_Array )" << "\n";
        output << "        {" << "\n";
        output << "            $i = $nNumber;" << "\n";
        output << "            $n = Floor( $nNumber / ( Count( $N_Array ) -1 ));" << "\n";
        output << "            $l = $nNumber - ( $n * ( Count( $N_Array ) -1 ));" << "\n";
        output << "" << "\n";
        output << "            Foreach( $N_Array as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                Unset( $R_Array[ $ACGT ][$N] );" << "\n";
        output << "                $R_Array[ $ACGT ][ $i != $l ? $i : 0 ] = $Value;" << "\n";
        output << "                $i -= $n;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== TempFile ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Order = [ 'GM', 'A', 'T', 'C', 'G' ];" << "\n";
        output << "" << "\n";
        output << "        $LTemp = Tempnam( '/tmp', 'L_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $MTemp = Tempnam( '/tmp', 'M_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $RTemp = Tempnam( '/tmp', 'R_'.$TSV_File.$Type_Name );" << "\n";
        output << "" << "\n";
        output << "        $lFtemp = Fopen( $LTemp, 'w' );" << "\n";
        output << "        $mFtemp = Fopen( $MTemp, 'w' );" << "\n";
        output << "        $rFtemp = Fopen( $RTemp, 'w' );" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $lFtemp, \"[\\n\" );" << "\n";
        output << "        Fwrite( $mFtemp, \"[\\n\" );" << "\n";
        output << "        Fwrite( $rFtemp, \"[\\n\" );" << "\n";
        output << "" << "\n";
        output << "        $Count1 = 0;" << "\n";
        output << "        For( $i = 0; $i < Count( $Order ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "" << "\n";
        output << "            if( !Array_Key_Exists( $Order[$i], $L_Array ))" << "\n";
        output << "                continue;" << "\n";
        output << "" << "\n";
        output << "            if( $Count1 != 0 ) Fwrite( $lFtemp, \",\\n\" );" << "\n";
        output << "            Fwrite( $lFtemp, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $L_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Count2 != 0 ) Fwrite( $lFtemp, ', ' );" << "\n";
        output << "                Fwrite( $lFtemp, '[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $lFtemp, \"]\\n    }\" );" << "\n";
        output << "            $Count1++;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Count1 = 0;" << "\n";
        output << "        For( $i = 0; $i < Count( $Order ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "" << "\n";
        output << "            if( !Array_Key_Exists( $Order[$i], $M_Array ))" << "\n";
        output << "                continue;" << "\n";
        output << "" << "\n";
        output << "            if( $Count1 != 0 ) Fwrite( $mFtemp, \",\\n\" );" << "\n";
        output << "            Fwrite( $mFtemp, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $M_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Count2 != 0 ) Fwrite( $mFtemp, ', ' );" << "\n";
        output << "                Fwrite( $mFtemp, '[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $mFtemp, \"]\\n    }\" );" << "\n";
        output << "            $Count1++;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Count1 = 0;" << "\n";
        output << "        For( $i = 0; $i < Count( $Order ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "" << "\n";
        output << "            if( !Array_Key_Exists( $Order[$i], $R_Array ))" << "\n";
        output << "                continue;" << "\n";
        output << "" << "\n";
        output << "            if( $Count1 != 0 ) Fwrite( $rFtemp, \",\\n\" );" << "\n";
        output << "            Fwrite( $rFtemp, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $R_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Count2 != 0 ) Fwrite( $rFtemp, ', ' );" << "\n";
        output << "                Fwrite( $rFtemp, '[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $rFtemp, \"]\\n    }\" );" << "\n";
        output << "            $Count1++;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $lFtemp, \"\\n]\\n\" );" << "\n";
        output << "        Fwrite( $mFtemp, \"\\n]\\n\" );" << "\n";
        output << "        Fwrite( $rFtemp, \"\\n]\\n\" );" << "\n";
        output << "" << "\n";
        output << "        Fclose( $lFtemp );" << "\n";
        output << "        Fclose( $mFtemp );" << "\n";
        output << "        Fclose( $rFtemp );" << "\n";
        output << "" << "\n";
        output << "#<!--================== SA_Plot ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Color_Array = $Ranking != 'Tailing％' ? \"['#FF0000','#088A08','#0000FF','#FFBF00']\" : \"['#000000','#FF0000','#088A08','#0000FF','#FFBF00']\";" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "            var margin = {top: 10, right: 20, bottom: 10, left: 20}" << "\n";
        output << "                width = ( svg_width - margin.left - margin.right ) / 5," << "\n";
        output << "                height = svg_height - margin.top  - margin.bottom;" << "\n";
        output << "" << "\n";
        output << "            var svgL = d3.select('body').append('div')" << "\n";
        output << "                .attr('id', 'svgL')" << "\n";
        output << "                .style('width', width * 2 - margin.right + 'px')" << "\n";
        output << "                .style('height', height + 'px')" << "\n";
        output << "                .style('display', 'inline-block' )" << "\n";
        output << "                .append('svg');" << "\n";
        output << "" << "\n";
        output << "            var svgM = d3.select('body').append('div')" << "\n";
        output << "                .attr('id', 'svgM')" << "\n";
        output << "                .style('width', width + margin.right + 5 + 'px')" << "\n";
        output << "                .style('height', height + 'px')" << "\n";
        output << "                .style('display', 'inline-block' )" << "\n";
        output << "                .append('svg');" << "\n";
        output << "" << "\n";
        output << "            var svgR = d3.select('body').append('div')" << "\n";
        output << "                .attr('id', 'svgR')" << "\n";
        output << "                .style('width', width * 2 - margin.right + 'px')" << "\n";
        output << "                .style('height', height + 'px')" << "\n";
        output << "                .style('display', 'inline-block' )" << "\n";
        output << "                .append('svg');" << "\n";
        output << "" << "\n";
        output << "            d3.json('$LTemp', function(data) {" << "\n";
        output << "                nv.addGraph(function() {" << "\n";
        output << "" << "\n";
        output << "                    var chart = nv.models.stackedAreaChart()" << "\n";
        output << "                        .x(function(d) { return d[0] })" << "\n";
        output << "                        .y(function(d) { return d[1] })" << "\n";
        output << "                        .color( $Color_Array )" << "\n";
        output << "                        .showControls(false)" << "\n";
        output << "                        .style('expand');" << "\n";
        output << "" << "\n";
        output << "                    d3.select('#svgL svg')" << "\n";
        output << "                      .datum(data)" << "\n";
        output << "                      .call(chart);" << "\n";
        output << "" << "\n";
        output << "                    nv.utils.windowResize(chart.update);" << "\n";
        output << "                    return chart;" << "\n";
        output << "                });" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            d3.json('$MTemp', function(data) {" << "\n";
        output << "                nv.addGraph(function() {" << "\n";
        output << "" << "\n";
        output << "                    var chart = nv.models.stackedAreaChart()" << "\n";
        output << "                        .x(function(d) { return d[0] })" << "\n";
        output << "                        .y(function(d) { return d[1] })" << "\n";
        output << "                        .color( $Color_Array )" << "\n";
        output << "                        .showControls(false)" << "\n";
        output << "                        .showYAxis(false)" << "\n";
        output << "                        .showXAxis(false)" << "\n";
        output << "                        .style('expand');" << "\n";
        output << "" << "\n";
        output << "                    d3.select('#svgM svg')" << "\n";
        output << "                      .datum(data)" << "\n";
        output << "                      .call(chart);" << "\n";
        output << "" << "\n";
        output << "                    nv.utils.windowResize(chart.update);" << "\n";
        output << "                    return chart;" << "\n";
        output << "                });" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            d3.json('$RTemp', function(data) {" << "\n";
        output << "                nv.addGraph(function() {" << "\n";
        output << "" << "\n";
        output << "                    var chart = nv.models.stackedAreaChart()" << "\n";
        output << "                        .margin({right: margin.right +12 })" << "\n";
        output << "                        .x(function(d) { return d[0] })" << "\n";
        output << "                        .y(function(d) { return d[1] })" << "\n";
        output << "                        .color( $Color_Array )" << "\n";
        output << "                        .rightAlignYAxis(true)" << "\n";
        output << "                        .showControls(false)" << "\n";
        output << "                        .xDomain([ $nNumber, 0 ])" << "\n";
        output << "                        .style('expand');" << "\n";
        output << "" << "\n";
        output << "                    d3.select('#svgR svg')" << "\n";
        output << "                      .datum(data)" << "\n";
        output << "                      .call(chart);" << "\n";
        output << "" << "\n";
        output << "                    nv.utils.windowResize(chart.update);" << "\n";
        output << "                    return chart;" << "\n";
        output << "                });" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=least>Least $Ranking</div>' );" << "\n";
        output << "            $( '#least' ).css({" << "\n";
        output << "                'position': 'absolute'," << "\n";
        output << "                'top': '35px'," << "\n";
        output << "                'left': '70px'," << "\n";
        output << "                'font-size': '16px'," << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=most>Most $Ranking</div>' );" << "\n";
        output << "            $( '#most' ).css({" << "\n";
        output << "                'position': 'absolute'," << "\n";
        output << "                'top': '35px'," << "\n";
        output << "                'left': width * 3 + 70 + 'px'," << "\n";
        output << "                'font-size': '16px'," << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=all>$Label</div>' );" << "\n";
        output << "            $( '#all' ).css({" << "\n";
        output << "                'position': 'absolute'," << "\n";
        output << "                'top': height - 15 + 'px'," << "\n";
        output << "                'left': width * 2 + width / $Label_Witdh + 'px'," << "\n";
        output << "                'font-size': '16px'," << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "        </script>\";" << "\n";
        output << "" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }

};

std::vector< std::map< std::string, SA_Type >> GeneTypeAnalyzerSA_Plot::anno_sa_table;

} // end of namespace algorithm
} // end of namespace ago
