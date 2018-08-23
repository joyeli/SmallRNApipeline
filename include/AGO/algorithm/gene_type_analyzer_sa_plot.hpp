#pragma once

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerSA_Plot
{
    struct SA_Type
    {
        std::vector< double > GMPM; // GMPM, GM, PM
        std::vector< double > Tail; // Atail, Ctail, Gtail, Other
        std::vector< double > tail; // aTail, cTail, gTail, Other is splite in to each Tail
        std::vector< std::map< char, double >> Ends; // 5'End-NT + Seed + 3'End-NT ( 0 ~ -4 ) + 3'End-4N-tail + 3'End-4N-notail ( 1 + 7 + 5 + 1 + 1 )
    
        SA_Type()
            : GMPM( 3, 0.0 )
            , Tail( 5, 0.0 )
            , tail( 4, 0.0 )
            , Ends( 15, std::map< char, double >() )
        {}
    };

  public:

    std::vector< std::map< std::string, SA_Type >> anno_sa_table;

    GeneTypeAnalyzerSA_Plot()
    {}

    static void make_sa_plot_table(
            const std::string& biotype,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< BedSampleType >& bed_samples,
            auto& genome_table,
            const std::size_t filter_ppm,
            std::vector< std::map< std::string, SA_Type >>& anno_sa_table,
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
            make_sa_counting_table( biotype, bed_samples[ smp ], anno_sa_table[ smp ], genome_table, filter_ppm, isSeed );
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
        char nt = ' ';
        std::string last4n_tail = "";
        std::string last4n_notail = "";

        sa_it->second.GMPM[0] += ppm;

        if( tail != 5 )
        {
            sa_it->second.GMPM[2] += ppm;
            sa_it->second.Tail[ tail ] += ppm;

            if( tail == 4 )
            {
                std::map< char, std::size_t > tail_count;

                for( std::size_t i = 0; i < tail_seq.length(); ++i )
                {
                    if( tail_count.find( tail_seq.at(i) ) == tail_count.end() )
                        tail_count[ tail_seq.at(i) ] = 0;
                    tail_count[ tail_seq.at(i) ] += 1;
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

            for( std::size_t i = 0; i < 5; ++i )
            {
                nt = sequence.at( sequence.length() -1 -i );

                if( md_map.find( sequence.length() -1 -i ) != md_map.end() ) nt = md_map[ sequence.length() -1 -i ];
                if( tc_set.find( sequence.length() -1 -i ) != tc_set.end() ) nt = 'T';

                if( sa_it->second.Ends[ 8 + i ].find( nt ) == sa_it->second.Ends[ 8 + i ].end() )
                    sa_it->second.Ends[ 8 + i ][ nt ] = 0;

                sa_it->second.Ends[ 8 + i ][ nt ] += ppm;
            }

            last4n_tail = sequence.substr( sequence.length() -1 -4, 4 );

            for( std::size_t i = 0; i < 4; ++i )
            {
                nt = last4n_tail.at(i);

                if( md_map.find( sequence.length() -1 -i ) != md_map.end() ) nt = md_map[ sequence.length() -1 -i ];
                if( tc_set.find( sequence.length() -1 -i ) != tc_set.end() ) nt = 'T';

                if( sa_it->second.Ends[13].find( nt ) == sa_it->second.Ends[13].end() )
                    sa_it->second.Ends[13][ nt ] = 0;

                sa_it->second.Ends[13][ nt ] += ppm;
            }
        }
        else
        {
            last4n_notail = sequence.substr( sequence.length() -1 -4, 4 );

            for( std::size_t i = 0; i < 4; ++i )
            {
                nt = last4n_notail.at(i);

                if( md_map.find( sequence.length() -1 -i ) != md_map.end() ) nt = md_map[ sequence.length() -1 -i ];
                if( tc_set.find( sequence.length() -1 -i ) != tc_set.end() ) nt = 'T';

                if( sa_it->second.Ends[14].find( nt ) == sa_it->second.Ends[14].end() )
                    sa_it->second.Ends[14][ nt ] = 0;

                sa_it->second.Ends[14][ nt ] += ppm;
            }

            sa_it->second.GMPM[1] += ppm;
        }

        for( std::size_t i = 0; i < 8; ++i )
        {
            nt = sequence.at(i);

            if( md_map.find(i) != md_map.end() ) nt = md_map[i];
            if( tc_set.find(i) != tc_set.end() ) nt = 'T';

            if( sa_it->second.Ends[i].find( nt ) == sa_it->second.Ends[i].end() )
                sa_it->second.Ends[i][ nt ] = 0;

            sa_it->second.Ends[i][ nt ] += ppm;
        }
    }

    static void make_sa_counting_table(
            const std::string& biotype,
            BedSampleType& bed_sample,
            std::map< std::string, SA_Type >& sa_table,
            auto& genome_table,
            const std::size_t filter_ppm,
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
            if( raw_bed.ppm_ < filter_ppm ) continue;
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

        for( std::size_t i = 0; i < 15; ++i ) sum.Ends[i]['S'] = 0.0;

        for( auto& sa : sa_table )
        {
            for( std::size_t i = 0; i < 5; ++i ) sum.Tail[i] += sa.second.Tail[i];
            for( std::size_t i = 0; i < 4; ++i ) sum.tail[i] += sa.second.tail[i];
            for( std::size_t i = 0; i < 15; ++i )
            for( auto & nt : sa.second.Ends[i] ) sum.Ends[i]['S'] += nt.second;
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
            for( std::size_t i = 0; i < 15; ++i )
            {
                // for( auto & nt : sa.second.Ends[i] ) nt.second = nt.second / sum.Ends[i]['S'];
                if( sa.second.Ends[i].find( 'A' ) == sa.second.Ends[i].end() ) sa.second.Ends[i][ 'A' ] = 0.0;
                if( sa.second.Ends[i].find( 'C' ) == sa.second.Ends[i].end() ) sa.second.Ends[i][ 'C' ] = 0.0;
                if( sa.second.Ends[i].find( 'G' ) == sa.second.Ends[i].end() ) sa.second.Ends[i][ 'G' ] = 0.0;
                if( sa.second.Ends[i].find( 'T' ) == sa.second.Ends[i].end() ) sa.second.Ends[i][ 'T' ] = 0.0;
            }
        }
    }

    static void output_sa_plot(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::map< std::string, SA_Type >& sa_table,
            const std::string& sample_name
            )
    {
        std::vector< std::ofstream > outputs( 17 );

        outputs[ 0].open( output_name + sample_name + "-tail.tsv" );
        outputs[ 1].open( output_name + sample_name + "-tail_with_other.tsv" );
        outputs[ 2].open( output_name + sample_name + "-nt_0.tsv" );
        outputs[ 3].open( output_name + sample_name + "-nt_1.tsv" );
        outputs[ 4].open( output_name + sample_name + "-nt_2.tsv" );
        outputs[ 5].open( output_name + sample_name + "-nt_3.tsv" );
        outputs[ 6].open( output_name + sample_name + "-nt_4.tsv" );
        outputs[ 7].open( output_name + sample_name + "-nt_5.tsv" );
        outputs[ 8].open( output_name + sample_name + "-nt_6.tsv" );
        outputs[ 9].open( output_name + sample_name + "-nt_7.tsv" );
        outputs[10].open( output_name + sample_name + "-nt_last_0.tsv" );
        outputs[11].open( output_name + sample_name + "-nt_last_-1.tsv" );
        outputs[12].open( output_name + sample_name + "-nt_last_-2.tsv" );
        outputs[13].open( output_name + sample_name + "-nt_last_-3.tsv" );
        outputs[14].open( output_name + sample_name + "-nt_last_-4.tsv" );
        outputs[15].open( output_name + sample_name + "-nt_last_4N-tail.tsv" );
        outputs[16].open( output_name + sample_name + "-nt_last_4N-notail.tsv" );

        for( std::size_t i = 0; i < outputs.size(); ++i )
        {
            outputs[i] << "Annotation\tGMPM\tGM\tPM\tA\tC\tG\tT";
            if( i == 1 ) outputs[i] << "\tO";
        }

        for( auto& anno : ano_len_idx.first )
        {
            if( sa_table.find( anno ) == sa_table.end() ) continue;

            for( auto& output : outputs ) output << "\n" << anno;

            for( std::size_t i = 0; i < 3; ++i )
                for( auto& output : outputs ) output << "\t" << sa_table[ anno ].GMPM[i];

            for( std::size_t i = 0; i < 4; ++i ) outputs[0] << "\t" << sa_table[ anno ].tail[i];
            for( std::size_t i = 0; i < 5; ++i ) outputs[1] << "\t" << sa_table[ anno ].Tail[i];

            for( std::size_t i = 0; i < 15; ++i )
                for( auto& nt : sa_table[ anno ].Ends[i] ) outputs[ i+2 ] << "\t" << nt.second;
        }

        for( auto& output : outputs ) output.close();
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
            for( std::size_t i = 0; i < 15; ++i )
            {
                std::cerr << ( sa.Ends[i].find( 'A' ) == sa.Ends[i].end() ? "\tAnotFound" : ( "\t" + std::to_string( sa.Ends[i][ 'A' ] )));
                std::cerr << ( sa.Ends[i].find( 'C' ) == sa.Ends[i].end() ? "\tCnotFound" : ( "\t" + std::to_string( sa.Ends[i][ 'C' ] )));
                std::cerr << ( sa.Ends[i].find( 'G' ) == sa.Ends[i].end() ? "\tGnotFound" : ( "\t" + std::to_string( sa.Ends[i][ 'G' ] )));
                std::cerr << ( sa.Ends[i].find( 'T' ) == sa.Ends[i].end() ? "\tTnotFound" : ( "\t" + std::to_string( sa.Ends[i][ 'T' ] )));
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
        output << "        $isGM = $_POST['isGM'];" << "\n";
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
        output << "        $Type_List = array('5’End','Seed0','Seed1','Seed2','Seed3','Seed4','Seed5','Seed6','3’End-4_Tailed','3’End-3_Tailed','3’End-2_Tailed','3’End-1_Tailed','3’End-0_Tailed','3’End4N_Tailed','3’End4N_NoTail','Tail','TailwithOther');" << "\n";
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
        output << "            case '5’End'          : $Type_Name = '-nt_0.tsv'              ; break;" << "\n";
        output << "            case 'Seed0'          : $Type_Name = '-nt_1.tsv'              ; break;" << "\n";
        output << "            case 'Seed1'          : $Type_Name = '-nt_2.tsv'              ; break;" << "\n";
        output << "            case 'Seed2'          : $Type_Name = '-nt_3.tsv'              ; break;" << "\n";
        output << "            case 'Seed3'          : $Type_Name = '-nt_4.tsv'              ; break;" << "\n";
        output << "            case 'Seed4'          : $Type_Name = '-nt_5.tsv'              ; break;" << "\n";
        output << "            case 'Seed5'          : $Type_Name = '-nt_6.tsv'              ; break;" << "\n";
        output << "            case 'Seed6'          : $Type_Name = '-nt_7.tsv'              ; break;" << "\n";
        output << "            case '3’End-4_Tailed' : $Type_Name = '-nt_last_-4.tsv'        ; break;" << "\n";
        output << "            case '3’End-3_Tailed' : $Type_Name = '-nt_last_-3.tsv'        ; break;" << "\n";
        output << "            case '3’End-2_Tailed' : $Type_Name = '-nt_last_-2.tsv'        ; break;" << "\n";
        output << "            case '3’End-1_Tailed' : $Type_Name = '-nt_last_-1.tsv'        ; break;" << "\n";
        output << "            case '3’End-0_Tailed' : $Type_Name = '-nt_last_0.tsv'         ; break;" << "\n";
        output << "            case '3’End4N_Tailed' : $Type_Name = '-nt_last_4N-tail.tsv'   ; break;" << "\n";
        output << "            case '3’End4N_NoTail' : $Type_Name = '-nt_last_4N-notail.tsv' ; break;" << "\n";
        output << "            case 'Tail'           : $Type_Name = '-tail.tsv'              ; break;" << "\n";
        output << "            case 'TailwithOther'  : $Type_Name = '-tail_with_other.tsv'   ; break;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "                    $Header[0] = $inFile_Line[2];" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 4; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                        $Header[ $i - 3 ] = $inFile_Line[$i];" << "\n";
        output << "" << "\n";
        output << "                    $isHeader = false;" << "\n";
        output << "                    continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Filter != '' && $Filter != 'FilterGMPM' && $inFile_Line[1] < $Filter  ) continue;" << "\n";
        output << "                if( $Ranking == 'Tailing％' && $inFile_Line[3] == 0 ) continue;" << "\n";
        output << "                if( $Ranking == 'GM' && $inFile_Line[2] == 0 ) continue;" << "\n";
        output << "                if( $Ranking == 'PM' && $inFile_Line[3] == 0 ) continue;" << "\n";
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
        output << "                $Ratio_Array[0] = $inFile_Line[2];" << "\n";
        output << "" << "\n";
        output << "                For( $i = 4; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    $Ratio_Array[ $i - 3 ] = $inFile_Line[$i];" << "\n";
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
        output << "" << "\n";
        output << "        if( $Type == '3’End4N_NoTail' )" << "\n";
        output << "            $Rank_List = array( 'GM' );" << "\n";
        output << "" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "            <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== isGM =====================-->" << "\n";
        output << "" << "\n";
        output << "        $GMbool = $isGM == 'True' ? true : false;" << "\n";
        output << "" << "\n";
        output << "        if( $Type != '5’End' && $Type != '3’End4N_NoTail' && Substr( $Type, 0, 1 ) != 'S' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<select name=isGM onchange=this.form.submit();>';" << "\n";
        output << "            echo '<option '; if($isGM=='') echo 'selected'; echo 'value= >isGM</option>';" << "\n";
        output << "" << "\n";
        output << "            $isGM_List = array( 'True', 'False' );" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $isGM_List ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$isGM_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $isGM == $isGM_List[$i] )" << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "                echo '>'.$isGM_List[$i].'</option>';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>" << "\n";
        output << "                <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "                <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "                <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "                <input type='hidden' name='Ranking' value='$Ranking' />" << "\n";
        output << "                <input type='hidden' name='nNumber' value='$nNumber' />" << "\n";
        output << "                <input type='hidden' name='isRatio' value='$isRatio' />" << "\n";
        output << "                <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "                <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "        }" << "\n";
        output << "        else $GMbool = false;" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "            <input type='hidden' name='isGM' value='$isGM' />" << "\n";
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
        output << "                    if( $isRatio == 'Yes' ) For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $L_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
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
        output << "                    if( $isRatio == 'Yes' ) For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $M_Array[ $Header[$j] ][0]  = $SumArray[$j];" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $M_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
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
        output << "                    if( $isRatio == 'Yes' ) For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $PM += $Temp_Array[$j][$k];" << "\n";
        output << "" << "\n";
        output << "                    For( $k = $GMbool ? 0 : 1; $k < Count( $Temp_Array[$j] ); ++$k )" << "\n";
        output << "                        $SumArray[$k] += $Temp_Array[$j][$k] / ( $isRatio == 'Yes' ? $PM : 1 );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / Count( $Temp_Array );" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $Sum += $SumArray[$j];" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $SumArray[$j] = $SumArray[$j] / $Sum;" << "\n";
        output << "                For( $j = $GMbool ? 0 : 1; $j < Count( $SumArray ); ++$j ) $R_Array[ $Header[$j] ][$N] = $SumArray[$j];" << "\n";
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
        output << "        $Order = $GMbool ? [ 'GM', 'A', 'T', 'C', 'G' ] : [ 'A', 'T', 'C', 'G' ];" << "\n";
        output << "" << "\n";
        output << "        $JsonL = Tempnam( '/tmp', 'Json_L_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $JsonM = Tempnam( '/tmp', 'Json_M_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $JsonR = Tempnam( '/tmp', 'Json_R_'.$TSV_File.$Type_Name );" << "\n";
        output << "" << "\n";
        output << "        $TextL = Tempnam( '/tmp', 'Text_L_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $TextM = Tempnam( '/tmp', 'Text_M_'.$TSV_File.$Type_Name );" << "\n";
        output << "        $TextR = Tempnam( '/tmp', 'Text_R_'.$TSV_File.$Type_Name );" << "\n";
        output << "" << "\n";
        output << "        $lJout = Fopen( $JsonL, 'w' );" << "\n";
        output << "        $mJout = Fopen( $JsonM, 'w' );" << "\n";
        output << "        $rJout = Fopen( $JsonR, 'w' );" << "\n";
        output << "" << "\n";
        output << "        $lTout = Fopen( $TextL, 'w' );" << "\n";
        output << "        $mTout = Fopen( $TextM, 'w' );" << "\n";
        output << "        $rTout = Fopen( $TextR, 'w' );" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $lJout, \"[\\n\" );" << "\n";
        output << "        Fwrite( $mJout, \"[\\n\" );" << "\n";
        output << "        Fwrite( $rJout, \"[\\n\" );" << "\n";
        output << "" << "\n";
        output << "        $Count1 = 0;" << "\n";
        output << "        For( $i = 0; $i < Count( $Order ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Count1 != 0 )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $lJout, \",\\n\" );" << "\n";
        output << "                Fwrite( $mJout, \",\\n\" );" << "\n";
        output << "                Fwrite( $rJout, \",\\n\" );" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $lJout, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "            Fwrite( $mJout, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "            Fwrite( $rJout, \"    {\\n        \\\"key\\\" : \\\"$Order[$i]\\\",\\n        \\\"values\\\" : [\" );" << "\n";
        output << "" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "            Foreach( $L_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $lJout, ( $Count2 == 0 ?  '' : ', ' ).'[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                if( $i != 0 ) $L_Array[ $Order[$i] ][$N] = $L_Array[ $Order[$i-1] ][$N] + $Value;" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "            Foreach( $M_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $mJout, ( $Count2 == 0 ?  '' : ', ' ).'[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                if( $i != 0 ) $M_Array[ $Order[$i] ][$N] = $M_Array[ $Order[$i-1] ][$N] + $Value;" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Count2 = 0;" << "\n";
        output << "            Foreach( $R_Array[ $Order[$i] ] as $N => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $rJout, ( $Count2 == 0 ?  '' : ', ' ).'[ '.$N.', '.$Value.' ]' );" << "\n";
        output << "                if( $i != 0 ) $R_Array[ $Order[$i] ][$N] = $R_Array[ $Order[$i-1] ][$N] + $Value;" << "\n";
        output << "                $Count2++;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $lJout, \"]\\n    }\" );" << "\n";
        output << "            Fwrite( $mJout, \"]\\n    }\" );" << "\n";
        output << "            Fwrite( $rJout, \"]\\n    }\" );" << "\n";
        output << "" << "\n";
        output << "            $Count1++;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $lJout, \"\\n]\\n\" );" << "\n";
        output << "        Fwrite( $mJout, \"\\n]\\n\" );" << "\n";
        output << "        Fwrite( $rJout, \"\\n]\\n\" );" << "\n";
        output << "" << "\n";
        output << "        Fclose( $lJout );" << "\n";
        output << "        Fclose( $mJout );" << "\n";
        output << "        Fclose( $rJout );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Order ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $i == 0 )" << "\n";
        output << "            {" << "\n";
        output << "                Foreach( $L_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $lTout, \"\\t\".$N );" << "\n";
        output << "                Foreach( $M_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $mTout, \"\\t\".$N );" << "\n";
        output << "                Foreach( $R_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $rTout, \"\\t\".$N );" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $lTout, \"\\n\".$Order[$i] );" << "\n";
        output << "            Fwrite( $mTout, \"\\n\".$Order[$i] );" << "\n";
        output << "            Fwrite( $rTout, \"\\n\".$Order[$i] );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $L_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $lTout, \"\\t\".$Value );" << "\n";
        output << "            Foreach( $M_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $mTout, \"\\t\".$Value );" << "\n";
        output << "            Foreach( $R_Array[ $Order[$i] ] as $N => $Value ) Fwrite( $rTout, \"\\t\".$Value );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fclose( $lTout );" << "\n";
        output << "        Fclose( $mTout );" << "\n";
        output << "        Fclose( $rTout );" << "\n";
        output << "" << "\n";
        output << "#<!--================== SA_Plot ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Color_Array = $GMbool ? \"['#000000','#FF0000','#088A08','#0000FF','#FFBF00']\" : \"['#FF0000','#088A08','#0000FF','#FFBF00']\";" << "\n";
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
        output << "            d3.json('$JsonL', function(data) {" << "\n";
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
        output << "            d3.json('$JsonM', function(data) {" << "\n";
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
        output << "            d3.json('$JsonR', function(data) {" << "\n";
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

} // end of namespace algorithm
} // end of namespace ago
