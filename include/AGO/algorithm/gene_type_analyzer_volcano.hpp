#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerVolcano
{

  public:

    GeneTypeAnalyzerVolcano()
    {}

    static void make_file_link( const std::string& p )
    {
        std::vector< std::string > list;
        boost::filesystem::path path( p + "../DiffBar/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
            if( file.path().filename().string().substr( 0, 20 ) == "LoadingDifferential_" )
                list.emplace_back( file.path().filename().string() );

        for( auto& file : list )
            if( !boost::filesystem::exists( p + file ))
                 boost::filesystem::create_symlink(( "../DiffBar/" + file ).c_str(), ( p + file ).c_str() );
    }

    static void output_volcano_visualization( const std::string& output_name, const std::string& biotype, const bool& isSeed )
    {
        make_file_link( output_name );
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $GMPM = $_POST['GMPM'];" << "\n";
        output << "        $isLog = $_POST['isLog'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $IsomiRs = $_POST['IsomiRs'];" << "\n";
        output << "        $Sample1 = $_POST['Sample1'];" << "\n";
        output << "        $Sample2 = $_POST['Sample2'];" << "\n";
        output << "        $Color_Low = $_POST['Color_Low'];" << "\n";
        output << "        $Color_Hight = $_POST['Color_Hight'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js></script>';" << "\n";
        output << "" << "\n";
        output << "        echo '<style type=\"text/css\">" << "\n";
        output << "                div[id=\"Dot\"] {" << "\n";
        output << "                    position: absolute; " << "\n";
        output << "                    text-align: left; " << "\n";
        output << "                    width: 180px;" << "\n";
        output << "                    height: 90px; " << "\n";
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
        output << "        $GMPM_List = array('GMPM', 'GM', 'PM', 'Atail', 'Ctail', 'Gtail', 'Ttail', 'Other');" << "\n";
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
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
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
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
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
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
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
            output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
            output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
            output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
            output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
            output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
            output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
            output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
            output << "            </form>\";" << "\n";
        }
        else
        {
            output << "        $IsomiRs == 'No';" << "\n";
        }

        output << "" << "\n";
        output << "#<!--================== isLog ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
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
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Colors =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        if( $Color_Low == '' ) $Color_Low = 'WhiteSmoke';" << "\n";
        output << "        if( $Color_Hight == '' ) $Color_Hight = 'Black';" << "\n";
        output << "        " << "\n";
        output << "        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        " << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == '' ) $Filter = 'FilterGMPM';" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value='.$Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample2' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Switch =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Sample1' value='$Sample2' />" << "\n";
        output << "            <input type='hidden' name='Sample2' value='$Sample1' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='submit' value='Switch Samples' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--===================== Read File ======================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = '';" << "\n";
        output << "" << "\n";
        output << "        if( File_Exists( './LoadingDifferential_'.$Sample1.'_'.$Sample2.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text' ))" << "\n";
        output << "            $TSV_File =  './LoadingDifferential_'.$Sample1.'_'.$Sample2.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text';" << "\n";
        output << "" << "\n";
        output << "        if( File_Exists( './LoadingDifferential_'.$Sample2.'_'.$Sample1.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text' ))" << "\n";
        output << "            $TSV_File =  './LoadingDifferential_'.$Sample2.'_'.$Sample1.( $IsomiRs == 'Yes' ? '-isomiRs' : '' ).'.text';" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $TSV_File.'_'.$GMPM.'_'.$Filter );" << "\n";
        output << "        $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $Ftemp, \"miRNA\\tLog2Fold\\tLogPvalue\\tGMPM\\n\" );" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $Column = 0;" << "\n";
        output << "" << "\n";
        output << "        $minPPM = 10000;" << "\n";
        output << "        $maxPPM = 0;" << "\n";
        output << "" << "\n";
        output << "        $minX = 0;" << "\n";
        output << "        $maxX = 0;" << "\n";
        output << "" << "\n";
        output << "        $minY = 0;" << "\n";
        output << "        $maxY = 0;" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "            if( $Filter != 'FilterGMPM' && $inFile_Line[1] < $Filter ) continue;" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                for( $i = 2; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    if( $inFile_Line[$i] == $Sample1.':'.$Sample2.':'.$GMPM ) $Column = $i;" << "\n";
        output << "" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( $isLog != '' && $inFile_Line[1] != 0 )" << "\n";
        output << "               $inFile_Line[1] = Log( $inFile_Line[1], $isLog );" << "\n";
        output << "" << "\n";
        output << "            if( $minPPM > $inFile_Line[1] ) $minPPM = $inFile_Line[1];" << "\n";
        output << "            if( $maxPPM < $inFile_Line[1] ) $maxPPM = $inFile_Line[1];" << "\n";
        output << "" << "\n";
        output << "            $Value = Explode( ':', $inFile_Line[ $Column ] );" << "\n";
        output << "" << "\n";
        output << "            $Xvalue =  Log( $Value[0], 2 );" << "\n";
        output << "            $Yvalue = -Log( $Value[1], 10 );" << "\n";
        output << "" << "\n";
        output << "            if( $maxX < Abs( $Xvalue )) $maxX = Abs( $Xvalue );" << "\n";
        output << "            if( $maxY < $Yvalue ) $maxY = $Yvalue;" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp," << "\n";
        output << "                $inFile_Line[0].\"\\t\"." << "\n";
        output << "                $Xvalue.\"\\t\"." << "\n";
        output << "                $Yvalue.\"\\t\"." << "\n";
        output << "                $inFile_Line[1].\"\\n\" );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $minX = -$maxX;" << "\n";
        output << "        Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "#<!--================== DotPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "            var color_map = d3.scale.linear()" << "\n";
        output << "                .domain([ $minPPM, $maxPPM ])" << "\n";
        output << "                .range([ '$Color_Low', '$Color_Hight' ]);" << "\n";
        output << "" << "\n";
        output << "            var margin = {top: 20, right: 120, bottom: 70, left: 40}" << "\n";
        output << "                width = svg_width - margin.left - margin.right," << "\n";
        output << "                height = svg_height - margin.top - margin.bottom;" << "\n";
        output << "" << "\n";
        output << "            var x = d3.scale.linear()" << "\n";
        output << "                .range([0, width]);" << "\n";
        output << "" << "\n";
        output << "            var y = d3.scale.linear()" << "\n";
        output << "                .range([height, 0]);" << "\n";
        output << "" << "\n";
        output << "            var color = d3.scale.category10();" << "\n";
        output << "" << "\n";
        output << "            var xAxis = d3.svg.axis()" << "\n";
        output << "                .scale(x)" << "\n";
        output << "                .orient('bottom');" << "\n";
        output << "" << "\n";
        output << "            var yAxis = d3.svg.axis()" << "\n";
        output << "                .scale(y)" << "\n";
        output << "                .orient('left');" << "\n";
        output << "" << "\n";
        output << "            var svg = d3.select('body').append('svg')" << "\n";
        output << "                .attr('width', width+70)" << "\n";
        output << "                .attr('height', height+60)" << "\n";
        output << "            .append('g')" << "\n";
        output << "                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n";
        output << "" << "\n";
        output << "            d3.tsv('$Temp', function(error, data) {" << "\n";
        output << "" << "\n";
        output << "                var xExtent = d3.extent(data, function(d) { return d.Log2Fold; });" << "\n";
        output << "                var yExtent = d3.extent(data, function(d) { return d.LogPvalue; });" << "\n";
        output << "" << "\n";
        output << "                xExtent[0] = $minX;" << "\n";
        output << "                yExtent[0] = $minY;" << "\n";
        output << "" << "\n";
        output << "                xExtent[1] = $maxX;" << "\n";
        output << "                yExtent[1] = $maxY;" << "\n";
        output << "" << "\n";
        output << "                x.domain( xExtent ).nice();" << "\n";
        output << "                y.domain( yExtent ).nice();" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'x axis')" << "\n";
        output << "                    .attr('transform', 'translate(0,' + height + ')')" << "\n";
        output << "                    .call(xAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('x', width)" << "\n";
        output << "                    .attr('y', -6)" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('Log2($Sample1/$Sample2)');" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'y axis')" << "\n";
        output << "                    .call(yAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('transform', 'rotate(-90)')" << "\n";
        output << "                    .attr('y', 6)" << "\n";
        output << "                    .attr('dy', '.71em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('-Log(Pvalue)');" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'xlineP')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x(1) )" << "\n";
        output << "                    .attr('y1', y(0) )" << "\n";
        output << "                    .attr('x2', x(1) )" << "\n";
        output << "                    .attr('y2', y($maxY) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'xlineN')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x(-1) )" << "\n";
        output << "                    .attr('y1', y(0) )" << "\n";
        output << "                    .attr('x2', x(-1) )" << "\n";
        output << "                    .attr('y2', y($maxY) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'yline1')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x($minX) )" << "\n";
        output << "                    .attr('y1', y(2) )" << "\n";
        output << "                    .attr('x2', x($maxX) )" << "\n";
        output << "                    .attr('y2', y(2) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "" << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'yline2')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x($minX) )" << "\n";
        output << "                    .attr('y1', y(1.3) )" << "\n";
        output << "                    .attr('x2', x($maxX) )" << "\n";
        output << "                    .attr('y2', y(1.3) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "" << "\n";
        output << "                var div = d3.select('body').append('div')" << "\n";
        output << "                    .attr('class', 'tooltip')" << "\n";
        output << "                    .attr('id', 'Dot')" << "\n";
        output << "                    .style('opacity', 0);" << "\n";
        output << "" << "\n";
        output << "                svg.selectAll('.dot')" << "\n";
        output << "                    .data(data)" << "\n";
        output << "                    .enter()" << "\n";

        if( !isSeed )
        {
            output << "                    .append('a')" << "\n";
            output << "                    .attr('xlink:href', function(d){ mirAnno = d.miRNA.split( '_' )[0]; return '../SqAlign/index.php?TSV_File=$Sample1.tsv&Annotation_Select=' + " << ( biotype == "miRNA" || biotype == "mirtron" || biotype == "miRNA_mirtron" ? "mirAnno.substring( 0, mirAnno.length -3 )" : "mirAnno" ) << "; })" << "\n";
            output << "                    .attr('target', '_blank')" << "\n";
        }
        else
        {
            output << "                    .append('a')" << "\n";
            output << "                    .attr('xlink:href', function(d){ mirAnno = d.miRNA.split( '_' )[0]; return '../SeedPie/index.php?TSV_File=All.tsv&Annotation_Select=' + mirAnno; })" << "\n";
            output << "                    .attr('target', '_blank')" << "\n";
        }

        output << "                    .append('circle')" << "\n";
        output << "                    .attr('class', 'dot')" << "\n";
        output << "                    .attr('r', 3)" << "\n";
        output << "                    .attr('cx', function(d) { return x(d.Log2Fold); })" << "\n";
        output << "                    .attr('cy', function(d) { return y(d.LogPvalue); })" << "\n";
        output << "                    .style('fill', function(d) { return color_map( d.GMPM ); })" << "\n";
        output << "                    .on('mouseover', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(200)" << "\n";
        output << "                            .style('opacity', .9);" << "\n";
        output << "" << "\n";
        output << "                        div.html('<table align=center ><tr><th colspan=2 >' + d.miRNA + '</th></tr><tr><th>TotalPPM</th><th>' +" << "\n";
        output << "                                 parseFloat( d.GMPM ).toFixed(4) + '</th></tr><tr><th>Log2(Fold)</th><th>' +" << "\n";
        output << "                                 parseFloat( d.Log2Fold ).toFixed(4) + '</th></tr><tr><th>-Log(Pvalue)</th><th>' +" << "\n";
        output << "                                 parseFloat( d.LogPvalue ).toFixed(2) + '</th></tr></table>')" << "\n";
        output << "                            .style('left', (d3.event.pageX) + 'px')" << "\n";
        output << "                            .style('top', (d3.event.pageY - 28) + 'px');" << "\n";
        output << "                    })" << "\n";
        output << "                    .on('mouseout', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(500)" << "\n";
        output << "                            .style('opacity', 0);" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    window.onresize = function(){" << "\n";
        output << "                    };" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=colorScale></div>' );" << "\n";
        output << "            $( '#colorScale' ).css({" << "\n";
        output << "                    'width': '7px'," << "\n";
        output << "                    'height': height + 'px'," << "\n";
        output << "                    'background': 'linear-gradient( to top, $Color_Low 0%, $Color_Hight 100% )'," << "\n";
        output << "                    'position': 'absolute'," << "\n";
        output << "                    'left': width + 55 + 'px'," << "\n";
        output << "                    'top': 50 + 'px'," << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "            for( var i = 100; i >= 0; i-=25 )" << "\n";
        output << "            {" << "\n";
        output << "                $( '#colorScale' ).append( '<div id=q' + i + '>-' + ( $minPPM + ( $maxPPM - $minPPM ) * ( i / 100 )).toFixed(2) + '</div>' );" << "\n";
        output << "                $( '#q' + i ).css({" << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': height * (( 100 - i )/100 ) - 10 + 'px'," << "\n";
        output << "                        'left': '7px'," << "\n";
        output << "                        });" << "\n";
        output << "            }" << "\n";
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
