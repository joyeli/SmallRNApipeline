#pragma once

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerSeedMap
{
  public:

    GeneTypeAnalyzerSeedMap()
    {}

    static double get_sum( std::map< std::string, double >& anno_ppm )
    {
        double sum = 0.0;
        for( auto& anno : anno_ppm ) sum += anno.second;
        return sum;
    }

    static void output_seed_mapping_table(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            std::vector< std::map< std::string, std::map< std::string, double >>>& seed_match_table
            )
    {
        std::ofstream output;
        std::set< std::string > anno_list;

        for( auto& seed_table : seed_match_table )
            for( auto& seed : seed_table )
                for( auto& anno : seed.second )
                    anno_list.emplace( anno.first );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            output.open( output_name + bed_samples[ smp ].first + ".tsv" );
            output << "Seed";

            for( auto& anno : anno_list )
                output << "\t" << anno;

            output << "\n";

            for( auto& seed : seed_match_table[ smp ] )
            {
                output << seed.first;

                for( auto& anno : anno_list )
                    output << "\t" << seed.second[ anno ];

                output << "\n";
            }
            output.close();
        }

        output.open( output_name + "All.tsv" );
        output << "Seed";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            output << "\t" << bed_samples[ smp ].first;

        output << "\n";

        for( auto& seed : seed_match_table[0] )
        {
            output << seed.first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                output << "\t" << get_sum( seed_match_table[ smp ][ seed.first ] );

            output << "\n";
        }
        output.close();
    }

    static void output_seedmap_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $isLog = $_POST['isLog'];" << "\n";
        output << "        $isBin = $_POST['isBin'];" << "\n";
        output << "        $Pivot = $_POST['Pivot'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $Color_Low = $_POST['Color_Low'];" << "\n";
        output << "        $Color_Pivot = $_POST['Color_Pivot'];" << "\n";
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
        output << "                    height: 70px; " << "\n";
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
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select Tsv</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size -1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "            if( $TSV_File == $TSV_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Pivot' value='$Pivot' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Pivot' value='$Color_Pivot' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
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
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='Pivot' value='$Pivot' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Pivot' value='$Color_Pivot' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Bin =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isBin onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isBin=='') echo 'selected'; echo 'value= >isBin</option>';" << "\n";
        output << "" << "\n";
        output << "        $isBin_List = array( 5, 10, 20, 40, 60, 80, 100 );" << "\n";
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
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='Pivot' value='$Pivot' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Color_Pivot' value='$Color_Pivot' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Colors =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $Pivotp = 0.5; " << "\n";
        output << "        $PivotP = 50; " << "\n";
        output << "" << "\n";
        output << "        if( $Color_Low == '' ) $Color_Low = 'White';" << "\n";
        output << "        if( $Color_Pivot == '' ) $Color_Pivot = 'WhiteSmoke';" << "\n";
        output << "        if( $Color_Hight == '' ) $Color_Hight = 'Black';" << "\n";
        output << "" << "\n";
        output << "        if( $Pivot == '' || $Pivot == \"Pivot%\" )" << "\n";
        output << "        {" << "\n";
        output << "            $Pivot = \"Pivot%\";" << "\n";
        output << "            $Pivotp = 0.5;" << "\n";
        output << "            $PivotP = 50;" << "\n";
        output << "        }" << "\n";
        output << "        else" << "\n";
        output << "        {" << "\n";
        output << "            $Pivotp = $Pivot / 100;" << "\n";
        output << "            $PivotP = $Pivot;" << "\n";
        output << "        }" << "\n";
        output << "        " << "\n";
        output << "        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Pivot size=8 value='.$Color_Pivot.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Pivot size=4 value='.$Pivot.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == '' ) $Filter = 'FilterGMPM';" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value='.$Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='isBin' value='$isBin' />" << "\n";
        output << "            <input type='hidden' name='isLog' value='$isLog' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Read Filter =====================-->" << "\n";
        output << "" << "\n";
        output << "        $Drop_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        if( $Filter != '' && $Filter != 'FilterGMPM' )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = new SplFileObject( 'All.tsv' );" << "\n";
        output << "            $isHeader = true;" << "\n";
        output << "" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "" << "\n";
        output << "                if( $inFile_Lines == '' ) continue;" << "\n";
        output << "                if( $isHeader )" << "\n";
        output << "                {" << "\n";
        output << "                    $isHeader = false;" << "\n";
        output << "                    continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "                $Sum = 0.0;" << "\n";
        output << "" << "\n";
        output << "                for( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    $Sum += Rtrim( $inFile_Line[$i] );" << "\n";
        output << "" << "\n";
        output << "                if( $Sum < ( Count( $inFile_Line ) -1 ) * $Filter )" << "\n";
        output << "                    Array_Push( $Drop_Array, $inFile_Line[0] );" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--===================== Read File ======================-->" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $TSV_File.'_'.$isLog );" << "\n";
        output << "        $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "" << "\n";
        output << "        $minPPM = 0;" << "\n";
        output << "        $maxPPM = 0;" << "\n";
        output << "" << "\n";
        output << "        $countX = 0;" << "\n";
        output << "        $countY = 0;" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                Fwrite( $Ftemp, \"X\\tY\\tSeed\\tPPM\\tAnno\\n\" );" << "\n";
        output << "" << "\n";
        output << "                for( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    $Anno_Array[$i] = Rtrim( $inFile_Line[$i] );" << "\n";
        output << "" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( In_Array( $inFile_Line[0], $Drop_Array ))" << "\n";
        output << "                continue;" << "\n";
        output << "" << "\n";
        output << "            $countX = 0;" << "\n";
        output << "            $countY++;" << "\n";
        output << "" << "\n";
        output << "            for( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $countX++;" << "\n";
        output << "                $Value = 0;" << "\n";
        output << "" << "\n";
        output << "                $topAnno = '';" << "\n";
        output << "                $topValue = 0;" << "\n";
        output << "" << "\n";
        output << "                for( $j = 0; $j < ( $isBin == '' ? 1 : $isBin ); ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $inFile_Line[ $i + $j ] = Rtrim( $inFile_Line[ $i + $j ]);" << "\n";
        output << "                    if( $topValue < $inFile_Line[ $i + $j ]) $topAnno = $Anno_Array[ $i + $j ];" << "\n";
        output << "                    $Value += $inFile_Line[ $i + $j ];" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $isBin != '' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Value = $Value / $isBin;" << "\n";
        output << "                    $i = $i + $isBin -1;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $isLog != '' && $Value != 0) $Value = Log( $Value, $isLog );" << "\n";
        output << "                if( $minPPM > $Value ) $minPPM = $Value;" << "\n";
        output << "                if( $maxPPM < $Value ) $maxPPM = $Value;" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, $countX.\"\\t\".$countY.\"\\t\".$inFile_Line[0].\"\\t\".$Value.\"\\t\".$topAnno.\"\\n\" );" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "#<!--================== DotPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "            var color_map = d3.scale.linear()" << "\n";
        output << "                .domain([ 0, $Pivotp, 1 ])" << "\n";
        output << "                .range([ '$Color_Low', '$Color_Pivot', '$Color_Hight' ]);" << "\n";
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
        output << "                var xExtent = d3.extent(data, function(d) { return d.X; });" << "\n";
        output << "                var yExtent = d3.extent(data, function(d) { return d.Y; });" << "\n";
        output << "" << "\n";
        output << "                xExtent[0] = 0;" << "\n";
        output << "                yExtent[0] = 0;" << "\n";
        output << "" << "\n";
        output << "                xExtent[1] = $countX;" << "\n";
        output << "                yExtent[1] = $countY;" << "\n";
        output << "" << "\n";
        output << "                x.domain( xExtent ).nice();" << "\n";
        output << "                y.domain( yExtent ).nice();" << "\n";
        output << "" << "\n";
        output << "                var div = d3.select('body').append('div')" << "\n";
        output << "                    .attr('id', 'Dot')" << "\n";
        output << "                    .style('opacity', 0);" << "\n";
        output << "" << "\n";
        output << "                svg.selectAll('.dot')" << "\n";
        output << "                    .data(data)" << "\n";
        output << "                    .enter()" << "\n";
        output << "                    .append('a')" << "\n";
        output << "                    .attr('xlink:href', function(d){ return '../SeedPie/index.php?TSV_File=$TSV_File&isLog=$isLog&Filter=$Filter&Annotation_Select=' + d.Seed })" << "\n";
        output << "                    .attr('target', '_blank')" << "\n";
        output << "                    .append('rect')" << "\n";
        output << "                    .attr('class', 'dot')" << "\n";
        output << "                    .attr('width', (x(2)-x(1)))" << "\n";
        output << "                    .attr('height',(y(1)-y(2)))" << "\n";
        output << "                    .attr('x', function(d) { return x(d.X-1); })" << "\n";
        output << "                    .attr('y', function(d) { return y(d.Y-1); })" << "\n";
        output << "                    .style('fill', function(d) { return color_map(d.PPM/$maxPPM); })" << "\n";
        output << "                    .on('mouseover', function(d) { if( d.PPM != 0 ) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(200)" << "\n";
        output << "                            .style('opacity', .9);" << "\n";
        output << "" << "\n";
        output << "                        div.html('<table align=center ><tr><th colspan=2 >' + d.Seed" << "\n";
        output << "                            + '</th></tr><tr><th>TopsAnno</th><th>' + d.Anno" << "\n";
        output << "                            + '</th></tr><tr><th>TotalPPM</th><th>' + parseFloat( d.PPM ).toFixed(4)" << "\n";
        output << "                            + '</th></tr></table>')" << "\n";
        output << "                            .style('left', (d3.event.pageX) + 'px')" << "\n";
        output << "                            .style('top', (d3.event.pageY - 28) + 'px');" << "\n";
        output << "                    }})" << "\n";
        output << "                    .on('mouseout', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(500)" << "\n";
        output << "                            .style('opacity', 0);" << "\n";
        output << "                    });" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=colorScale></div>' );" << "\n";
        output << "            $( '#colorScale' ).css({" << "\n";
        output << "                    'width': '7px'," << "\n";
        output << "                    'height': height + 'px'," << "\n";
        output << "                    'background': 'linear-gradient( to top, $Color_Low 0%, $Color_Pivot $PivotP%, $Color_Hight 100% )'," << "\n";
        output << "                    'position': 'absolute'," << "\n";
        output << "                    'left': width + 55 + 'px'," << "\n";
        output << "                    'top': 50 + 'px'," << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "            for( var i = 100; i >= 0; i-=25 )" << "\n";
        output << "            {" << "\n";
        output << "                $( '#colorScale' ).append( '<div id=q' + i + '>' + ( $minPPM + ( $maxPPM - $minPPM ) * ( i / 100 )).toFixed(2) + '</div>' );" << "\n";
        output << "                $( '#q' + i ).css({" << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': height * (( 100 - i )/100 ) - 10 + 'px'," << "\n";
        output << "                        'left': '7px'," << "\n";
        output << "                        });" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( $Pivotp != 0.5 )" << "\n";
        output << "            {" << "\n";
        output << "                $( '#colorScale' ).append( '<div id=colorTop >' + ( $minPPM + ( $maxPPM - $minPPM ) * $Pivotp ).toFixed(2) + '</div>' );" << "\n";
        output << "                $( '#colorTop' ).css({" << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': height * ( 1 - $Pivotp ) - 10 + 'px'," << "\n";
        output << "                        'color': 'red'," << "\n";
        output << "                        'left': '7px'," << "\n";
        output << "                        });" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=textY>Seeds</div>' );" << "\n";
        output << "            $( '#textY' ).css({" << "\n";
        output << "                'position': 'absolute'," << "\n";
        output << "                'top': '100px'," << "\n";
        output << "                'left': '-10px'," << "\n";
        output << "                'font-size': '25px'," << "\n";
        output << "                'transform': 'rotatez(-90deg)'" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "            $( 'body' ).append( '<div id=textX>Annotations</div>' );" << "\n";
        output << "            $( '#textX' ).css({" << "\n";
        output << "                'position': 'absolute'," << "\n";
        output << "                'top': '30px'," << "\n";
        output << "                'left': '50px'," << "\n";
        output << "                'font-size': '25px'," << "\n";
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
