#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerValplot
{

  public:

    GeneTypeAnalyzerValplot()
    {}

    static void output_valplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            const std::string& token
            )
    {
        std::ofstream output( output_name + token + ".tsv" );
        output << "Annotation";

        double gm = 0.0;
        double pm = 0.0;

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : ano_len_idx.first )
        {
            output << "\n" << anno << std::setprecision( 0 ) << std::fixed
                << ( anno_mark[0].find( anno ) != anno_mark[0].end() ? anno_mark[0][ anno ] : "" );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                gm = 0.0;
                pm = 0.0;

                if( token != "PM" )
                {
                    if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                                gm += anno_table_tail[ smp ][5][ anno ][ len ];
                        }
                }

                if( token != "GM" )
                {
                    for( std::size_t i = 0; i < 5; i++ )
                    {
                        if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                            for( auto& len : ano_len_idx.second )
                            {
                                if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                    pm += anno_table_tail[ smp ][i][ anno ][ len ];
                            }
                    }
                }

                output << "\t" << ( token == "GMPM" ? gm + pm : ( token == "GM" ? gm : ( token == "PM" ? pm : (( gm + pm ) < 1 ? 0 : ( pm * 100 / ( gm + pm ))))));
            }
        }

        output << "\n";
        output.close();
    }

    static void output_valplot_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n"; 
        output << "<html>" << "\n"; 
        output << "    <meta charset='utf-8'>" << "\n"; 
        output << "    <body>" << "\n"; 
        output << "" << "\n"; 
        output << "    <? " << "\n"; 
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n"; 
        output << "        $FGMPM = $_POST['FGMPM'];" << "\n"; 
        output << "        $isLog = $_POST['isLog'];" << "\n"; 
        output << "        $Filter = $_POST['Filter'];" << "\n"; 
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n"; 
        output << "        $FilterMin = $_POST['FilterMin'];" << "\n"; 
        output << "        $FilterMax = $_POST['FilterMax'];" << "\n"; 
        output << "        $isAbundant = $_POST['isAbundant'];" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n"; 
        output << "        echo '<script src=http://code.jquery.com/jquery-1.7.min.js ></script>';" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/jquery.event.drag-2.0.min.js ></script>';" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/slick.core.js ></script>';" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/slick.grid.js ></script>';" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/slick.dataview.js ></script>';" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/d3.parcoords.js ></script>';" << "\n"; 
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/divgrid.js ></script>';" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<link rel=stylesheet type=text/css href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/slick.grid.css />';" << "\n"; 
        output << "        echo '<link rel=stylesheet type=text/css href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/jquery-ui-1.8.16.custom.css />';" << "\n"; 
        output << "        echo '<link rel=stylesheet type=text/css href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/slickgrid/examples.css />';" << "\n"; 
        output << "        echo '<link rel=stylesheet type=text/css href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/d3.parcoords.css />';" << "\n"; 
        output << "        echo '<link rel=stylesheet type=text/css href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/style.css />';" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<style type=\"text/css\">" << "\n"; 
        output << "                html, body {" << "\n"; 
        output << "                    height: 99%;" << "\n"; 
        output << "                    width: 99%;" << "\n"; 
        output << "                }" << "\n"; 
        output << "                " << "\n"; 
        output << "                #grid {" << "\n"; 
        output << "                    position: fixed;" << "\n"; 
        output << "                    width: 99%;" << "\n"; 
        output << "                    bottom: 0;" << "\n"; 
        output << "                    height: 300px;" << "\n"; 
        output << "                }" << "\n"; 
        output << "            </style>';" << "\n"; 
        output << "" << "\n"; 
        output << "#<!--================== TSV File ====================-->" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n"; 
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n"; 
        output << "" << "\n"; 
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n"; 
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n"; 
        output << "        $List_Size = Count( $TSV_List );" << "\n"; 
        output << "" << "\n"; 
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n"; 
        output << "        {" << "\n"; 
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n"; 
        output << "" << "\n"; 
        output << "            if( $TSV_File == $TSV_List[$i] ) " << "\n"; 
        output << "                echo 'selected ';" << "\n"; 
        output << "" << "\n"; 
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n"; 
        output << "        }" << "\n"; 
        output << "        echo \"</select>" << "\n"; 
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n"; 
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n"; 
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n"; 
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n"; 
        output << "            </form>\";" << "\n"; 
        output << "" << "\n"; 
        output << "#<!--================== is_Abundant ====================-->" << "\n"; 
        output << "                " << "\n"; 
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n"; 
        output << "        echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n"; 
        output << "        " << "\n"; 
        output << "        $isAbundant_List = array('MostAbundant', 'AllAnnotation');" << "\n"; 
        output << "        $isAbundant_Size = Count( $isAbundant_List );" << "\n"; 
        output << "        if( $isAbundant == '' )" << "\n"; 
        output << "            $isAbundant = 'MostAbundant';" << "\n"; 
        output << "        " << "\n"; 
        output << "        For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n"; 
        output << "        {" << "\n"; 
        output << "            echo '<option value='.$isAbundant_List[$i].' ';" << "\n"; 
        output << "        " << "\n"; 
        output << "            if( $isAbundant == $isAbundant_List[$i] )" << "\n"; 
        output << "                echo 'selected ';" << "\n"; 
        output << "        " << "\n"; 
        output << "            echo '>' . $isAbundant_List[$i] . '</option>';" << "\n"; 
        output << "        }" << "\n"; 
        output << "" << "\n"; 
        output << "        echo \"</select>" << "\n"; 
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n"; 
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n"; 
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n"; 
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n"; 
        output << "            </form>\";" << "\n"; 
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
        output << "            <input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n"; 
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n"; 
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n"; 
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n"; 
        output << "            </form>\";" << "\n"; 
        output << "" << "\n"; 
        output << "#<!--================== Filter Min & Max ====================-->" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n"; 
        output << "        echo '<input type=text name=FilterMin size=5 value=';" << "\n"; 
        output << "" << "\n"; 
        output << "        if( $FilterMin=='' )" << "\n"; 
        output << "            echo 'FilterMin';" << "\n"; 
        output << "        else" << "\n"; 
        output << "            echo $FilterMin;" << "\n"; 
        output << "" << "\n"; 
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n"; 
        output << "        echo '<input type=text name=FilterMax size=5 value=';" << "\n"; 
        output << "" << "\n"; 
        output << "        if( $FilterMax=='' )" << "\n"; 
        output << "            echo 'FilterMax';" << "\n"; 
        output << "        else" << "\n"; 
        output << "            echo $FilterMax;" << "\n"; 
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
        output << "#<!--================== Filter GMPM ====================-->" << "\n"; 
        output << "" << "\n"; 
        output << "        echo '<select name=FGMPM >';" << "\n"; 
        output << "        echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';" << "\n"; 
        output << "" << "\n"; 
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n"; 
        output << "        {" << "\n"; 
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n"; 
        output << "" << "\n"; 
        output << "            if( $FGMPM == $TSV_List[$i] ) " << "\n"; 
        output << "                echo 'selected ';" << "\n"; 
        output << "" << "\n"; 
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n"; 
        output << "        }" << "\n"; 
        output << "" << "\n"; 
        output << "        echo \"" << "\n"; 
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n"; 
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n"; 
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n"; 
        output << "            <input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n"; 
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n"; 
        output << "            <input type='submit' value='Submit' /> " << "\n"; 
        output << "            </form><br/>\";" << "\n"; 
        output << "" << "\n"; 
        output << "#<!--================== Read File ====================-->" << "\n"; 
        output << "" << "\n"; 
        output << "        if( $TSV_File != '' )" << "\n"; 
        output << "        {" << "\n"; 
        output << "        " << "\n"; 
        output << "            if( $FGMPM != '' && $Filter != '' )" << "\n"; 
        output << "            {" << "\n"; 
        output << "                $FFile = $FGMPM;" << "\n"; 
        output << "                $FinFile = File_get_contents( $FFile );" << "\n"; 
        output << "                $FinFile_Lines = Explode( \"\\n\", $FinFile );" << "\n"; 
        output << "                $FArray = Array();" << "\n"; 
        output << "            }" << "\n"; 
        output << "" << "\n"; 
        output << "            $inFile = File_get_contents( $TSV_File );" << "\n"; 
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n"; 
        output << "" << "\n"; 
        output << "            $Temp = Tempnam( '/tmp', $TSV_File );" << "\n"; 
        output << "            $Ftemp = Fopen( $Temp, 'w' );" << "\n"; 
        output << "" << "\n"; 
        output << "            $Sample_Nu = Count( Explode( \"\\t\", $inFile_Lines[0] ))-1;" << "\n"; 
        output << "            $Value_Max = 0;" << "\n"; 
        output << "            $Value_Min = 0;" << "\n"; 
        output << "" << "\n"; 
        output << "            For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n"; 
        output << "            {" << "\n"; 
        output << "                if( $FGMPM != '' && $Filter != '' && $i != 0 )" << "\n"; 
        output << "                {" << "\n"; 
        output << "                    $Check_Filter = Array();" << "\n"; 
        output << "                    $FinFile_Line = Explode( \"\\t\", $FinFile_Lines[$i] );" << "\n"; 
        output << "" << "\n"; 
        output << "                    For( $j = 1; $j < Count( $FinFile_Line ); ++$j )" << "\n"; 
        output << "                    {" << "\n"; 
        output << "                        $FValue = $FinFile_Line[$j];" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( $isLog != '' && $FValue != 0 )" << "\n"; 
        output << "                            $FValue = ( Log($FinFile_Line[$j]) / Log($isLog) );" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( $FValue <= $Filter )" << "\n"; 
        output << "                            Array_Push( $Check_Filter, $FValue );" << "\n"; 
        output << "                    }" << "\n"; 
        output << "" << "\n"; 
        output << "                    if( Count( $Check_Filter ) == Count( $FinFile_Line )-1 )" << "\n"; 
        output << "                        Continue;" << "\n"; 
        output << "                }" << "\n"; 
        output << "" << "\n"; 
        output << "" << "\n"; 
        output << "                $Check_Filter = 'keep';" << "\n"; 
        output << "" << "\n"; 
        output << "                if( $i != 0 )" << "\n"; 
        output << "                {" << "\n"; 
        output << "                    $Check_Min = Array();" << "\n"; 
        output << "                    $Check_Max = Array();" << "\n"; 
        output << "" << "\n"; 
        output << "                    $inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n"; 
        output << "" << "\n"; 
        output << "                    For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n"; 
        output << "                    {" << "\n"; 
        output << "                        $Value = $inFile_Line[$j];" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( $isLog != '' && $Value != 0 )" << "\n"; 
        output << "                            $Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( $FilterMin != 'FilterMin' && $FilterMin != '' )" << "\n"; 
        output << "                            if( $Value <= $FilterMin )" << "\n"; 
        output << "                                Array_Push( $Check_Min, $Value );" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( $FilterMax != 'FilterMax' && $FilterMax != '' )" << "\n"; 
        output << "                            if( $Value >= $FilterMax )" << "\n"; 
        output << "                                Array_Push( $Check_Max, $Value );" << "\n"; 
        output << "                    }" << "\n"; 
        output << "" << "\n"; 
        output << "                    if( Count( $Check_Min ) == Count( $inFile_Line )-1 || Count( $Check_Max ) == Count( $inFile_Line )-1 )" << "\n"; 
        output << "                        $Check_Filter = 'notkeep';" << "\n"; 
        output << "" << "\n"; 
        output << "" << "\n"; 
        output << "                    if( $isAbundant == 'MostAbundant' )" << "\n"; 
        output << "                    {" << "\n"; 
        output << "                        $miRNA = Explode( '*', $inFile_Line[0] );" << "\n"; 
        output << "" << "\n"; 
        output << "                        if( Count( $miRNA ) != 2 )" << "\n"; 
        output << "                            $Check_Filter = 'notkeep';" << "\n"; 
        output << "                    }" << "\n"; 
        output << "" << "\n"; 
        output << "                    if( $Check_Filter == 'keep' )" << "\n"; 
        output << "                    {" << "\n"; 
        output << "                        For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n"; 
        output << "                        {" << "\n"; 
        output << "                            $Value = $inFile_Line[$j];" << "\n"; 
        output << "" << "\n"; 
        output << "                            if( $isLog != '' && $Value != 0 )" << "\n"; 
        output << "                                $Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n"; 
        output << "" << "\n"; 
        output << "                            if( $Value_Max <= $Value )" << "\n"; 
        output << "                                $Value_Max =  $Value;" << "\n"; 
        output << "" << "\n"; 
        output << "                            if( $Value_Min >= $Value )" << "\n"; 
        output << "                                $Value_Min =  $Value;" << "\n"; 
        output << "                        }" << "\n"; 
        output << "                    }" << "\n"; 
        output << "                }" << "\n"; 
        output << "" << "\n"; 
        output << "                if( $Check_Filter == 'keep' )" << "\n"; 
        output << "                {" << "\n"; 
        output << "                    if( $isLog != '' && $i != 0)" << "\n"; 
        output << "                    {" << "\n"; 
        output << "                        $inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n"; 
        output << "                        Fwrite( $Ftemp, $inFile_Line[0] );" << "\n"; 
        output << "" << "\n"; 
        output << "                        For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n"; 
        output << "                        {" << "\n"; 
        output << "                            if( $inFile_Line[$j] != 0 )" << "\n"; 
        output << "                            {    " << "\n"; 
        output << "                                $Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n"; 
        output << "                                Fwrite( $Ftemp, \"\\t\".$Value );" << "\n"; 
        output << "                            }" << "\n"; 
        output << "                            else" << "\n"; 
        output << "                                Fwrite( $Ftemp, \"\\t\".'0' );" << "\n"; 
        output << "                        }" << "\n"; 
        output << "" << "\n"; 
        output << "                        Fwrite( $Ftemp, \"\\n\" );" << "\n"; 
        output << "                    }" << "\n"; 
        output << "                    else" << "\n"; 
        output << "                        Fwrite( $Ftemp, $inFile_Lines[$i].\"\\n\" );" << "\n"; 
        output << "                }" << "\n"; 
        output << "            }" << "\n"; 
        output << "" << "\n"; 
        output << "            Fwrite( $Ftemp, 'MinValue' );" << "\n"; 
        output << "            For( $i = 0; $i < $Sample_Nu; ++$i )" << "\n"; 
        output << "                Fwrite( $Ftemp, \"\\t\".$Value_Min );" << "\n"; 
        output << "            Fwrite( $Ftemp, \"\\n\" );" << "\n"; 
        output << "" << "\n"; 
        output << "            Fwrite( $Ftemp, 'MaxValue' );" << "\n"; 
        output << "            For( $i = 0; $i < $Sample_Nu; ++$i )" << "\n"; 
        output << "                Fwrite( $Ftemp, \"\\t\".$Value_Max );" << "\n"; 
        output << "            Fwrite( $Ftemp, \"\\n\" );" << "\n"; 
        output << "" << "\n"; 
        output << "            fclose( $Ftemp );" << "\n"; 
        output << "" << "\n"; 
        output << "#<!--================== ValPlot ====================-->" << "\n"; 
        output << "    " << "\n"; 
        output << "            echo \"<div id='example' class='parcoords' style='height:240px;'></div>" << "\n"; 
        output << "                <div id='grid'></div>" << "\n"; 
        output << "                <script id='brushing'>" << "\n"; 
        output << "                    var parcoords = d3.parcoords()('#example')" << "\n"; 
        output << "                            .alpha(0.4)" << "\n"; 
        output << "                            .mode('queue') // progressive rendering" << "\n"; 
        output << "                            .height(d3.max([document.body.clientHeight-326, 220]))" << "\n"; 
        output << "                            .margin({" << "\n"; 
        output << "                                top: 36," << "\n"; 
        output << "                                left: 0," << "\n"; 
        output << "                                right: 0," << "\n"; 
        output << "                                bottom: 16" << "\n"; 
        output << "                            });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                    // load csv file and create the chart" << "\n"; 
        output << "                    d3.tsv('$Temp', function(data) {" << "\n"; 
        output << "                        // slickgrid needs each data element to have an id" << "\n"; 
        output << "                        data.forEach(function(d,i) { d.id = d.id || i; });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        parcoords" << "\n"; 
        output << "                            .data(data)" << "\n"; 
        output << "                            .hideAxis(['Annotation'])" << "\n"; 
        output << "                            .hideAxis(['id'])" << "\n"; 
        output << "                            .render()" << "\n"; 
        output << "                            .reorderable()" << "\n"; 
        output << "                            .brushMode('1D-axes');" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // setting up grid" << "\n"; 
        output << "                        var column_keys = d3.keys(data[0]);" << "\n"; 
        output << "                        var columns = column_keys.map(function(key,i) {" << "\n"; 
        output << "                            return {" << "\n"; 
        output << "                                id: key," << "\n"; 
        output << "                                name: key," << "\n"; 
        output << "                                field: key," << "\n"; 
        output << "                                sortable: true" << "\n"; 
        output << "                            }" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        var options = {" << "\n"; 
        output << "                            enableCellNavigation: true," << "\n"; 
        output << "                            enableColumnReorder: false," << "\n"; 
        output << "                            multiColumnSort: false" << "\n"; 
        output << "                        };" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        var dataView = new Slick.Data.DataView();" << "\n"; 
        output << "                        var grid = new Slick.Grid('#grid', dataView, columns, options);" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // wire up model events to drive the grid" << "\n"; 
        output << "                        dataView.onRowCountChanged.subscribe(function (e, args) {" << "\n"; 
        output << "                            grid.updateRowCount();" << "\n"; 
        output << "                            grid.render();" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        dataView.onRowsChanged.subscribe(function (e, args) {" << "\n"; 
        output << "                            grid.invalidateRows(args.rows);" << "\n"; 
        output << "                            grid.render();" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // column sorting" << "\n"; 
        output << "                        var sortcol = column_keys[0];" << "\n"; 
        output << "                        var sortdir = 1;" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        function comparer(a, b) {" << "\n"; 
        output << "                            var x = a[sortcol], y = b[sortcol];" << "\n"; 
        output << "                            return (x == y ? 0 : (x > y ? 1 : -1));" << "\n"; 
        output << "                        }" << "\n"; 
        output << "                        " << "\n"; 
        output << "                        // click header to sort grid column" << "\n"; 
        output << "                        grid.onSort.subscribe(function (e, args) {" << "\n"; 
        output << "                            sortdir = args.sortAsc ? 1 : -1;" << "\n"; 
        output << "                            sortcol = args.sortCol.field;" << "\n"; 
        output << "                    " << "\n"; 
        output << "                            if ($.browser.msie && $.browser.version <= 8) {" << "\n"; 
        output << "                                dataView.fastSort(sortcol, args.sortAsc);" << "\n"; 
        output << "                            } else {" << "\n"; 
        output << "                                dataView.sort(comparer, args.sortAsc);" << "\n"; 
        output << "                            }" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // highlight row in chart" << "\n"; 
        output << "                        grid.onMouseEnter.subscribe(function(e,args) {" << "\n"; 
        output << "                            var i = grid.getCellFromEvent(e).row;" << "\n"; 
        output << "                            var d = parcoords.brushed() || data;" << "\n"; 
        output << "                            parcoords.highlight([d[i]]);" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                        grid.onMouseLeave.subscribe(function(e,args) {" << "\n"; 
        output << "                            parcoords.unhighlight();" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // fill grid with data" << "\n"; 
        output << "                        gridUpdate(data);" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        // update grid on brush" << "\n"; 
        output << "                        parcoords.on('brush', function(d) {" << "\n"; 
        output << "                            gridUpdate(d);" << "\n"; 
        output << "                        });" << "\n"; 
        output << "                    " << "\n"; 
        output << "                        function gridUpdate(data) {" << "\n"; 
        output << "                            dataView.beginUpdate();" << "\n"; 
        output << "                            dataView.setItems(data);" << "\n"; 
        output << "                            dataView.endUpdate();" << "\n"; 
        output << "                        };" << "\n"; 
        output << "                    });" << "\n"; 
        output << "                </script>\";" << "\n"; 
        output << "        }" << "\n"; 
        output << "    ?>" << "\n"; 
        output << "    </body>" << "\n"; 
        output << "</html>" << "\n"; 

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
