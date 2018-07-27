#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDotplot
{

  public:

    GeneTypeAnalyzerDotplot()
    {}

    static void output_dotplot_isomirs(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::ofstream output( output_name + sample_name + "-isomiRs.tsv" );
        output << sample_name << "\tGMPM\tGM\tPM\tTailing_Ratio\tA_Tail\tC_Tail\tG_Tail\tT_Tail\tOther_Tail\n";

        double gm = 0.0;
        double pm = 0.0;
        double tail_a = 0.0;
        double tail_c = 0.0;
        double tail_g = 0.0;
        double tail_t = 0.0;
        double tail_o = 0.0;

        for( auto& anno : ano_len_idx.first )
        {
            gm = 0.0;
            pm = 0.0;
            tail_a = 0.0;
            tail_c = 0.0;
            tail_g = 0.0;
            tail_t = 0.0;
            tail_o = 0.0;

            if( anno_table_tail[5].find( anno ) != anno_table_tail[5].end() )
                for( auto& len : anno_table_tail[5][ anno ] ) gm += len.second;

            for( std::size_t i = 0; i < 5; i++ )
            {
                if( anno_table_tail[i].find( anno ) == anno_table_tail[i].end() ) continue;
                for( auto& len : anno_table_tail[i][ anno ] ){ pm += len.second;
                switch(i)
                {
                    case 0 : tail_a += len.second; break;
                    case 1 : tail_c += len.second; break;
                    case 2 : tail_g += len.second; break;
                    case 3 : tail_t += len.second; break;
                    case 4 : tail_o += len.second; break;
                }}
            }

            output
                << std::setprecision( 0 )
                << std::fixed
                << anno << ( anno_mark.find( anno ) != anno_mark.end() ? anno_mark[ anno ] : "" ) << "\t"
                << gm + pm << "\t"
                << gm << "\t"
                << pm << "\t"
                << (( gm + pm ) < 1 ? 0 :( pm * 100 /( gm + pm ))) << "\t"
                << tail_a << "\t"
                << tail_c << "\t"
                << tail_g << "\t"
                << tail_t << "\t"
                << tail_o << "\n"
                ;
        }

        output.close();
    }

    static void output_dotplot(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::string anno_temp;
        std::vector< std::string > split;

        std::set< std::string > anno_mark_set;
        std::map< std::string, std::vector< double >> anno_map;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));

            if( anno_mark.find( anno ) != anno_mark.end() )
            {
                if( anno_mark[ anno ] == "!" || anno_mark[ anno ] == "*!" )
                    anno_mark_set.emplace( split[0] );
            }
        }

        std::ofstream output( output_name + sample_name + ".tsv" );
        output << sample_name << "\tGMPM\tGM\tPM\tTailing_Ratio\tA_Tail\tC_Tail\tG_Tail\tT_Tail\tOther_Tail\n";

        double gm = 0.0;
        double pm = 0.0;
        double tail_a = 0.0;
        double tail_c = 0.0;
        double tail_g = 0.0;
        double tail_t = 0.0;
        double tail_o = 0.0;

        for( auto& anno : ano_len_idx.first )
        {
            gm = 0.0;
            pm = 0.0;
            tail_a = 0.0;
            tail_c = 0.0;
            tail_g = 0.0;
            tail_t = 0.0;
            tail_o = 0.0;

            if( anno_table_tail[5].find( anno ) != anno_table_tail[5].end() )
                for( auto& len : anno_table_tail[5][ anno ] ) gm += len.second;

            for( std::size_t i = 0; i < 5; i++ )
            {
                if( anno_table_tail[i].find( anno ) == anno_table_tail[i].end() ) continue;
                for( auto& len : anno_table_tail[i][ anno ] ){ pm += len.second;
                switch(i)
                {
                    case 0 : tail_a += len.second; break;
                    case 1 : tail_c += len.second; break;
                    case 2 : tail_g += len.second; break;
                    case 3 : tail_t += len.second; break;
                    case 4 : tail_o += len.second; break;
                }}
            }

            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            anno_temp = split[0] + ( anno_mark_set.find( anno ) == anno_mark_set.end() ? "" : "!" );

            if( anno_map.find( anno_temp ) == anno_map.end() )
                anno_map[ anno_temp ] = std::vector< double >( 7, 0.0 );

            anno_map[ anno_temp ][0] += gm;
            anno_map[ anno_temp ][1] += pm;
            anno_map[ anno_temp ][2] += tail_a;
            anno_map[ anno_temp ][3] += tail_c;
            anno_map[ anno_temp ][4] += tail_g;
            anno_map[ anno_temp ][5] += tail_t;
            anno_map[ anno_temp ][6] += tail_o;
        }

        for( auto& anno : anno_map )
        {
            gm     = anno.second[0];
            pm     = anno.second[1];
            tail_a = anno.second[2];
            tail_c = anno.second[3];
            tail_g = anno.second[4];
            tail_t = anno.second[5];
            tail_o = anno.second[6];

            output
                << std::setprecision( 0 )
                << std::fixed
                << anno.first << "\t"
                << gm + pm << "\t"
                << gm << "\t"
                << pm << "\t"
                << (( gm + pm ) < 1 ? 0 :( pm * 100 /( gm + pm ))) << "\t"
                << tail_a << "\t"
                << tail_c << "\t"
                << tail_g << "\t"
                << tail_t << "\t"
                << tail_o << "\n"
                ;
        }

        output.close();
    }

    static void output_dotplot_visualization( const std::string& output_name, const std::string& biotype, const bool& isSeed )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $GMPM = $_POST['GMPM'];" << "\n";
        output << "        $FGMPM = $_POST['FGMPM'];" << "\n";
        output << "        $isLog = $_POST['isLog'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $IsomiRs = $_POST['IsomiRs'];" << "\n";
        output << "        $ForceMin = $_POST['ForceMin'];" << "\n";
        output << "        $ForceMax = $_POST['ForceMax'];" << "\n";
        output << "        $TSV_File1 = $_POST['TSV_File1'];" << "\n";
        output << "        $TSV_File2 = $_POST['TSV_File2'];" << "\n";
        output << "        $Color_Low = $_POST['Color_Low'];" << "\n";
        output << "        $Color_Hight = $_POST['Color_Hight'];" << "\n";
        output << "" << "\n";
        output << "        $isAbundant = $IsomiRs == 'No' ? 'AllmiRNA' : $_POST['isAbundant'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js></script>';" << "\n";
        output << "" << "\n";
        output << "        echo '<style type=\"text/css\">" << "\n";
        output << "                div[id=\"Dot\"] {" << "\n";
        output << "                    position: absolute; " << "\n";
        output << "                    text-align: left; " << "\n";
        output << "                    width: 180px;" << "\n";
        output << "                    height: 80px; " << "\n";
        output << "                    padding: 8px; " << "\n";
        output << "                    background: lightsteelblue; " << "\n";
        output << "                    font: 12px sans-serif;" << "\n";
        output << "                    border: 0px;" << "\n";
        output << "                    border-radius: 8px; " << "\n";
        output << "                    pointer-events: none; " << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                th {" << "\n";
        output << "                    border: 1px solid black;" << "\n";
        output << "                } " << "\n";
        output << "" << "\n";
        output << "                .axis path, .axis line {" << "\n";
        output << "                    fill: none;" << "\n";
        output << "                    stroke: #000;" << "\n";
        output << "                    shape-rendering: crispEdges;" << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                #grid {" << "\n";
        output << "                    position: fixed;" << "\n";
        output << "                    width: 100%;" << "\n";
        output << "                    bottom: 0;" << "\n";
        output << "                    height: 300px;" << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                .slick-row:hover {" << "\n";
        output << "                    font-weight: bold;" << "\n";
        output << "                    color: #069;" << "\n";
        output << "                }" << "\n";
        output << "            </style>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=GMPM onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($GMPM=='') echo 'selected'; echo '>GM or PM</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPM_List = array('GMPM', 'GM', 'PM', 'Tailing_Ratio', 'A_Tail', 'C_Tail', 'G_Tail', 'T_Tail', 'Other_Tail');" << "\n";
        output << "" << "\n";
        output << "        $GMPM_Size = Count( $GMPM_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $GMPM_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$GMPM_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $GMPM == $GMPM_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $GMPM_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n";
        output << "            <input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n";
        output << "            <input type='hidden' name='TSV_File1' value='$TSV_File1' />" << "\n";
        output << "            <input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";

        if( !isSeed && biotype != "BioType/" )
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
            output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
            output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n";
            output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
            output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
            output << "            <input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n";
            output << "            <input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n";
            output << "            <input type='hidden' name='TSV_File1' value='$TSV_File1' />" << "\n";
            output << "            <input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n";
            output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
            output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
            output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
            output << "            </form>\";" << "\n";
        }
        else
        {
            output << "        $IsomiRs == 'No';" << "\n";
        }

        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "        $TSV_List_Temp = array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $TSV_Temp = Explode( '-isomiRs', $TSV_List[$i] );" << "\n";
        output << "" << "\n";
        output << "            if( $IsomiRs == 'Yes' )" << "\n";
        output << "            {" << "\n";
        output << "                if( Count( $TSV_Temp ) > 1 )" << "\n";
        output << "                    Array_Push( $TSV_List_Temp, $TSV_Temp[0] );" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "            {" << "\n";
        output << "                if( Count( $TSV_Temp ) == 1 )" << "\n";
        output << "                    Array_Push( $TSV_List_Temp, Substr( $TSV_Temp[0], 0, Strlen( $TSV_Temp[0] ) - 4 ));" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $TSV_List = $TSV_List_Temp;" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File1 onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File1=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File1 == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo '</select>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File2 onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File2=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File2 == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";

        if( biotype != "BioType/" )
        {
            output << "        $TSV_File_1 = $TSV_File1.( $IsomiRs == 'Yes' ? '-isomiRs.tsv' : '.tsv' );" << "\n";
            output << "        $TSV_File_2 = $TSV_File2.( $IsomiRs == 'Yes' ? '-isomiRs.tsv' : '.tsv' );" << "\n";
        }
        else
        {
            output << "        $TSV_File_1 = $TSV_File1.'-isomiRs.tsv';" << "\n";
            output << "        $TSV_File_2 = $TSV_File2.'-isomiRs.tsv';" << "\n";
        }

        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n";
        output << "            <input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";

        if( !isSeed && biotype != "BioType/" )
        {
            output << "#<!--================== is_Abundant ====================-->" << "\n";
            output << "" << "\n";
            output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
            output << "        echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n";
            output << "" << "\n";
            output << "        $isAbundant_List = array('MostAbundant', 'AllmiRNA');" << "\n";
            output << "        $isAbundant_Size = Count( $isAbundant_List );" << "\n";
            output << "" << "\n";
            output << "        if( $isAbundant == '' )" << "\n";
            output << "            $isAbundant = 'MostAbundant';" << "\n";
            output << "" << "\n";
            output << "        For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n";
            output << "        {" << "\n";
            output << "            echo '<option value='.$isAbundant_List[$i].' ';" << "\n";
            output << "" << "\n";
            output << "            if( $isAbundant == $isAbundant_List[$i] )" << "\n";
            output << "                echo 'selected ';" << "\n";
            output << "" << "\n";
            output << "            echo '>' . $isAbundant_List[$i] . '</option>';" << "\n";
            output << "        }" << "\n";
            output << "" << "\n";
            output << "        echo \"</select>" << "\n";
            output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
            output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n";
            output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
            output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
            output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
            output << "            <input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n";
            output << "            <input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n";
            output << "            <input type='hidden' name='TSV_File1' value='$TSV_File1' />" << "\n";
            output << "            <input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n";
            output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
            output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
            output << "            </form>\";" << "\n";
        }
        else
        {
            output << "        $isAbundant = 'AllmiRNA';" << "\n";
        }

        output << "" << "\n";
        output << "#<!--================== isLog ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=isLog onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isLog=='') echo 'selected'; echo 'value= >isLog</option>';" << "\n";
        output << "" << "\n";
        output << "        $isLog_List = array(2, 4, 6, 8, 10);" << "\n";
        output << "        $isLog_Size = Count( $isLog_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $isLog_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$isLog_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isLog == $isLog_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$isLog_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n";
        output << "            <input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n";
        output << "            <input type='hidden' name='TSV_File1' value='$TSV_File1' />" << "\n";
        output << "            <input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== ForceMin&Max ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=ForceMin size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $ForceMin=='' )" << "\n";
        output << "            echo 'ForceMin';" << "\n";
        output << "        else" << "\n";
        output << "            echo $ForceMin;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        echo '<input type=text name=ForceMax size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $ForceMax=='' )" << "\n";
        output << "            echo 'ForceMax';" << "\n";
        output << "        else" << "\n";
        output << "            echo $ForceMax;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Filter size=4 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Filter=='' )" << "\n";
        output << "            echo 'Filter';" << "\n";
        output << "        else" << "\n";
        output << "            echo $Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Colors =====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Color_Low == '' ) $Color_Low = 'WhiteSmoke';" << "\n";
        output << "        if( $Color_Hight == '' ) $Color_Hight = 'Black';" << "\n";
        output << "        " << "\n";
        output << "        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        " << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=FGMPM >';" << "\n";
        output << "        echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';" << "\n";
        output << "" << "\n";
        output << "        $FGMPM_List = array('GMPM', 'GM', 'PM', 'Tailing_Ratio', 'A_Tail', 'C_Tail', 'G_Tail', 'T_Tail', 'Other_Tail');" << "\n";
        output << "        $FGMPM_Size = Count( $FGMPM_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $FGMPM_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$FGMPM_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $FGMPM == $FGMPM_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $FGMPM_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='TSV_File1' value='$TSV_File1' />" << "\n";
        output << "            <input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Index & Header ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Files = array( $TSV_File_1, $TSV_File_2 );" << "\n";
        output << "        $Index = array();" << "\n";
        output << "        $Column = 0;" << "\n";
        output << "        $FColumn = 0;" << "\n";
        output << "" << "\n";
        output << "        For( $l = 0; $l < Count( $Files ); ++$l )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = File_get_contents( $Files[$l] );" << "\n";
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n";
        output << "" << "\n";
        output << "                if( $i == 0 )" << "\n";
        output << "                {" << "\n";
        output << "                    For( $j = 0; $j < Count( $inFile_Line ); ++$j )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( $inFile_Line[$j] == $GMPM )" << "\n";
        output << "                            $Column = $j;" << "\n";
        output << "" << "\n";
        output << "                        if( $inFile_Line[$j] == $FGMPM )" << "\n";
        output << "                            $FColumn = $j;" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "                else" << "\n";
        output << "                {" << "\n";
        output << "                    $miRNA_Seed = Explode( '*', $inFile_Line[0] );" << "\n";
        output << "" << "\n";
        output << "                    if( $isAbundant == 'MostAbundant' && $IsomiRs == 'Yes' && Count( $miRNA_Seed ) == 1 )" << "\n";
        output << "                        continue;" << "\n";
        output << "" << "\n";
        output << "                    Array_Push( $Index, $miRNA_Seed[0] );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Sort( $Index, SORT_STRING );" << "\n";
        output << "        $uIndex = Array_Unique( $Index, SORT_STRING );" << "\n";
        output << "" << "\n";
        output << "#<!--================== Read File & Log ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Anno_Value1 = array();" << "\n";
        output << "        $Anno_Value2 = array();" << "\n";
        output << "        $Filter_Array1 = array();" << "\n";
        output << "        $Filter_Array2 = array();" << "\n";
        output << "        $Sample_Name = array();" << "\n";
        output << "        $FSample_Name = array();" << "\n";
        output << "" << "\n";
        output << "        For( $l = 0; $l < Count( $Files ); ++$l )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = File_get_contents( $Files[$l] );" << "\n";
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n";
        output << "" << "\n";
        output << "                if( $i == 0 )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $l == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Anno_Value1['Anno'] = $inFile_Line[0];" << "\n";
        output << "                        $Filter_Array1['Anno']=$inFile_Line[0].'_F'.$FGMPM;" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $Anno_Value2['Anno'] = $inFile_Line[0];" << "\n";
        output << "                        $Filter_Array2['Anno']=$inFile_Line[0].'_F'.$FGMPM;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    Array_Push( $Sample_Name, $inFile_Line[0] );" << "\n";
        output << "                    Array_Push( $FSample_Name, $inFile_Line[0].'_F'.$FGMPM );" << "\n";
        output << "                }" << "\n";
        output << "                else" << "\n";
        output << "                {" << "\n";
        output << "                    $miRNA_Seed = Explode( '*', $inFile_Line[0] );" << "\n";
        output << "                    $Value = $inFile_Line[$Column];" << "\n";
        output << "                    $FValue = $inFile_Line[$FColumn];" << "\n";
        output << "" << "\n";
        output << "                    if( $isLog != '' && $Value != 0 )" << "\n";
        output << "                        $Value = ( Log($inFile_Line[$Column]) / Log($isLog) );" << "\n";
        output << "" << "\n";
        output << "                    if( $isLog != '' && $FValue != 0 )" << "\n";
        output << "                        $FValue = ( Log($inFile_Line[$FColumn]) / Log($isLog) );" << "\n";
        output << "" << "\n";
        output << "                    if( $l == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Anno_Value1[$miRNA_Seed[0]] = Round($Value,2);" << "\n";
        output << "                        $Filter_Array1[$miRNA_Seed[0]]=Round($FValue,2);" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $Anno_Value2[$miRNA_Seed[0]] = Round($Value,2);" << "\n";
        output << "                        $Filter_Array2[$miRNA_Seed[0]]=Round($FValue,2);" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter & Temp ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $GMPM.'_'.$Sample_Name[0].'_'.$Sample_Name[1].'_'.$isAbundant.'_'.$isLog.'_'.$Filter.'_'.$FGMPM );" << "\n";
        output << "        $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "        $MaxAxis = 0;" << "\n";
        output << "" << "\n";
        output << "        $minPPM = 1000000;" << "\n";
        output << "        $maxPPM = 0;" << "\n";
        output << "" << "\n";
        output << "        $Filter = $isLog != '' && $isLog != 'isLog' ? ( Log( $Filter ) / Log( $isLog )) : $Filter;" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $Ftemp, 'miRNA'.\"\\t\"." << "\n";
        output << "                $Anno_Value1['Anno'].\"\\t\"." << "\n";
        output << "                $Anno_Value2['Anno'].\"\\t\"." << "\n";
        output << "                $Filter_Array1['Anno'].\"\\t\"." << "\n";
        output << "                $Filter_Array2['Anno'].\"\\n\" );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Index ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $uIndex[$i] != '' && $uIndex[$i] != 'Anno' )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Anno_Value1[$uIndex[$i]] == '' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Anno_Value1[$uIndex[$i]]  = 0;" << "\n";
        output << "                    $Filter_Array1[$uIndex[$i]]= 0;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Anno_Value2[$uIndex[$i]] == '' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Anno_Value2[$uIndex[$i]]  = 0;" << "\n";
        output << "                    $Filter_Array2[$uIndex[$i]]= 0;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Filter != 'Filter' && $FGMPM != '' )" << "\n";
        output << "                    if( $Filter_Array1[$uIndex[$i]] < $Filter && $Filter_Array2[$uIndex[$i]] < $Filter )" << "\n";
        output << "                        continue;" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, $uIndex[$i].\"\\t\"." << "\n";
        output << "                                $Anno_Value1[$uIndex[$i]].\"\\t\"." << "\n";
        output << "                                $Anno_Value2[$uIndex[$i]].\"\\t\"." << "\n";
        output << "                                $Filter_Array1[$uIndex[$i]].\"\\t\"." << "\n";
        output << "                                $Filter_Array2[$uIndex[$i]].\"\\n\" );" << "\n";
        output << "" << "\n";
        output << "                if( $MaxAxis < $Anno_Value1[$uIndex[$i]] ) $MaxAxis = $Anno_Value1[$uIndex[$i]];" << "\n";
        output << "                if( $MaxAxis < $Anno_Value2[$uIndex[$i]] ) $MaxAxis = $Anno_Value2[$uIndex[$i]];" << "\n";
        output << "" << "\n";
        output << "                if( $maxPPM < $Filter_Array1[$uIndex[$i]] ) $maxPPM = $Filter_Array1[$uIndex[$i]];" << "\n";
        output << "                if( $maxPPM < $Filter_Array2[$uIndex[$i]] ) $maxPPM = $Filter_Array2[$uIndex[$i]];" << "\n";
        output << "                if( $minPPM > $Filter_Array1[$uIndex[$i]] ) $minPPM = $Filter_Array1[$uIndex[$i]];" << "\n";
        output << "                if( $minPPM > $Filter_Array2[$uIndex[$i]] ) $minPPM = $Filter_Array2[$uIndex[$i]];" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "        if( $FGMPM == '' )" << "\n";
        output << "            $FGMPM = 'Filter';" << "\n";
        output << "" << "\n";
        output << "        if( $ForceMin == 'ForceMin' || $ForceMin == '' )" << "\n";
        output << "            $ForceMin = 0;" << "\n";
        output << "" << "\n";
        output << "        if( $ForceMax == 'ForceMax' || $ForceMax == '' )" << "\n";
        output << "            $ForceMax = $MaxAxis;" << "\n";
        output << "" << "\n";
        output << "#<!--================== DotPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "        " << "\n";
        output << "            var margin = {top: 20, right: 70, bottom: 70, left: 70}" << "\n";
        output << "                width = svg_width - margin.left - margin.right," << "\n";
        output << "                height = svg_height - margin.top - margin.bottom;" << "\n";
        output << "        " << "\n";
        output << "            var color_map = d3.scale.linear()" << "\n";
        output << "                .domain([ $minPPM, $maxPPM ])" << "\n";
        output << "                .range([ '$Color_Low', '$Color_Hight' ]);" << "\n";
        output << "        " << "\n";
        output << "            var x = d3.scale.linear()" << "\n";
        output << "                .range([0, width]);" << "\n";
        output << "        " << "\n";
        output << "            var y = d3.scale.linear()" << "\n";
        output << "                .range([height, 0]);" << "\n";
        output << "        " << "\n";
        output << "            var color = d3.scale.category10();" << "\n";
        output << "        " << "\n";
        output << "            var xAxis = d3.svg.axis()" << "\n";
        output << "                .scale(x)" << "\n";
        output << "                .orient('bottom');" << "\n";
        output << "        " << "\n";
        output << "            var yAxis = d3.svg.axis()" << "\n";
        output << "                .scale(y)" << "\n";
        output << "                .orient('left');" << "\n";
        output << "        " << "\n";
        output << "            var svg = d3.select('body').append('svg')" << "\n";
        output << "                .attr('width', width+100)" << "\n";
        output << "                .attr('height', height+40)" << "\n";
        output << "            .append('g')" << "\n";
        output << "                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n";
        output << "        " << "\n";
        output << "            d3.tsv('$Temp', function(error, data) {" << "\n";
        output << "        " << "\n";
        output << "                var xExtent = d3.extent(data, function(d) { return d.$Sample_Name[1]; });" << "\n";
        output << "                var yExtent = d3.extent(data, function(d) { return d.$Sample_Name[0]; });" << "\n";
        output << "        " << "\n";
        output << "                xExtent[0] = $ForceMin;" << "\n";
        output << "                yExtent[0] = $ForceMin;" << "\n";
        output << "        " << "\n";
        output << "                xExtent[1] = $ForceMax;" << "\n";
        output << "                yExtent[1] = $ForceMax;" << "\n";
        output << "        " << "\n";
        output << "                x.domain( xExtent ).nice();" << "\n";
        output << "                y.domain( yExtent ).nice();" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'x axis')" << "\n";
        output << "                    .attr('transform', 'translate(0,' + height + ')')" << "\n";
        output << "                    .call(xAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('x', width)" << "\n";
        output << "                    .attr('y', -6)" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('$Sample_Name[1]');" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'y axis')" << "\n";
        output << "                    .call(yAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('transform', 'rotate(-90)')" << "\n";
        output << "                    .attr('y', 6)" << "\n";
        output << "                    .attr('dy', '.71em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('$Sample_Name[0]');" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'xyline')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x($ForceMin) )" << "\n";
        output << "                    .attr('y1', y($ForceMin) )" << "\n";
        output << "                    .attr('x2', x($ForceMax) )" << "\n";
        output << "                    .attr('y2', y($ForceMax) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "        " << "\n";
        output << "                var div = d3.select('body').append('div')" << "\n";
        output << "                    .attr('class', 'tooltip')" << "\n";
        output << "                    .attr('id', 'Dot')" << "\n";
        output << "                    .style('opacity', 0);" << "\n";
        output << "        " << "\n";
        output << "                svg.selectAll('.dot')" << "\n";
        output << "                    .data(data)" << "\n";
        output << "                    .enter()" << "\n";

        if( !isSeed && biotype != "BioType/" )
        {
            output << "                    .append('a')" << "\n";
        	output << "                    .attr('xlink:href', function(d){ mirAnno = d.miRNA.split( '_' )[0]; return '../SqAlign/index.php?TSV_File=$TSV_File1.tsv&Annotation_Select=' + " << ( biotype == "miRNA" || biotype == "mirtron" || biotype == "miRNA_mirtron" ? "mirAnno.substring( 0, mirAnno.length -3 )" : "mirAnno" ) << "; })" << "\n";
            output << "                    .attr('target', '_blank')" << "\n";
        }

        output << "                    .append('circle')" << "\n";
        output << "                    .attr('class', 'dot')" << "\n";
        output << "                    .attr('r', 5)" << "\n";
        output << "                    .attr('cx', function(d) { return x(d.$Sample_Name[1]); })" << "\n";
        output << "                    .attr('cy', function(d) { return y(d.$Sample_Name[0]); })" << "\n";
        output << "                    .style('fill', function(d) { return '$FGMPM' == 'Filter' ? color(d.species) : color_map( d.$Sample_Name[0]_F$FGMPM > d.$Sample_Name[1]_F$FGMPM ? d.$Sample_Name[0]_F$FGMPM : d.$Sample_Name[1]_F$FGMPM ); })" << "\n";
        output << "                    .on('mouseover', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(200)" << "\n";
        output << "                            .style('opacity', .9);" << "\n";
        output << "        " << "\n";
        output << "                        div.html('<table align=center ><tr><th colspan=3 >' + d.miRNA +" << "\n";
        output << "                                 '</th></tr><tr><th>Sample</th><th>$GMPM</th><th>$FGMPM</th></tr><tr><th>$Sample_Name[0]</th><th>' +" << "\n";
        output << "                                 d.$Sample_Name[0] + '</th><th>' + d.$Sample_Name[0]_F$FGMPM + '</th></tr><tr><th>$Sample_Name[1]</th><th>' +" << "\n";
        output << "                                 d.$Sample_Name[1] + '</th><th>' + d.$Sample_Name[1]_F$FGMPM + '</th></tr></table>')" << "\n";
        output << "                            .style('left', (d3.event.pageX) + 'px')" << "\n";
        output << "                            .style('top', (d3.event.pageY - 28) + 'px');" << "\n";
        output << "                    })" << "\n";
        output << "                    .on('mouseout', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(500)" << "\n";
        output << "                            .style('opacity', 0);" << "\n";
        output << "                    });" << "\n";
        output << "        " << "\n";
        output << "                    window.onresize = function(){" << "\n";
        output << "                    };" << "\n";
        output << "            });" << "\n";
        output << "        </script>\";" << "\n";
        output << "        " << "\n";
        output << "        " << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
