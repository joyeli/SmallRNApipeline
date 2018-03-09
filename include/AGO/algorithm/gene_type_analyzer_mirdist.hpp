#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerMirdist
{

  public:

    GeneTypeAnalyzerMirdist()
    {}

    static void output_mirdist(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >> anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            const std::string& token
            )
    {
        std::ofstream output;
        std::vector< std::pair< std::string, double >> total_sorted;

        double gm = 0.0;
        double pm = 0.0;

        for( auto& anno : ano_len_idx.first )
        {
            gm = 0.0;
            pm = 0.0;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                    for( auto& len : ano_len_idx.second )
                    {
                        if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                            gm += anno_table_tail[ smp ][5][ anno ][ len ];
                    }

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

            total_sorted.emplace_back( anno, gm + pm );
        }

        std::sort( total_sorted.begin(), total_sorted.end(),
        []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
        {
            if( a.second == b.second )
                return a.first > b.first;
            else
                return a.second > b.second;
        });

        for( auto& len : ano_len_idx.second )
        {
            output.open( output_name + token + "_" + std::to_string( len ) + ".tsv" );
            output << "Annotation";

            for( auto& smp  : bed_samples ) output << "\t" << smp.first;
            for( auto& anno : total_sorted )
            {
                output << "\n" << anno.first << std::setprecision( 0 ) << std::fixed;

                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                {
                    gm = 0.0;
                    pm = 0.0;

                    if( token != "PM" || token == "Tailing" )
                    {
                        if( anno_table_tail[ smp ][5].find( anno.first ) != anno_table_tail[ smp ][5].end() )
                            if( anno_table_tail[ smp ][5][ anno.first ].find( len ) != anno_table_tail[ smp ][5][ anno.first ].end() )
                                gm += anno_table_tail[ smp ][5][ anno.first ][ len ];
                    }

                    if( token != "GM" || token == "Tailing" )
                    {
                        for( std::size_t i = 0; i < 5; i++ )
                        {
                            if( anno_table_tail[ smp ][i].find( anno.first ) != anno_table_tail[ smp ][i].end() )
                                if( anno_table_tail[ smp ][i][ anno.first ].find( len ) != anno_table_tail[ smp ][i][ anno.first ].end() )
                                    pm += anno_table_tail[ smp ][i][ anno.first ][ len ];
                        }
                    }

                    output << "\t" << ( token == "GMPM" ? gm + pm : ( token == "GM" ? gm : ( token == "PM" ? pm : (( gm + pm ) < 1 ? 0 : ( pm * 100 / ( gm + pm ))))));
                }
            }

            output << "\n";
            output.close();
        }
    }

    static void output_mirdist_visualization(
            const std::string& output_name
            )
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
        output << "        $MaxHight = $_POST['MaxHight'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $Top_miRNA = $_POST['Top_miRNA'];" << "\n";
        output << "        $mirDistType = $_POST['mirDistType'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/d3.min.js></script>';" << "\n";
        output << "        echo '<link href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/svg0331.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=GMPM onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($GMPM=='') echo 'selected'; echo '>GM or PM</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPM_List = array('GMPM', 'GM', 'PM', 'Tailing');" << "\n";
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
        output << "            <input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='mirDistType' value='$mirDistType' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "        echo '<style>" << "\n";
        output << "                .x.axis path {" << "\n";
        output << "                display: none;" << "\n";
        output << "            }" << "\n";
        output << "            </style>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Top miRNA Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Top_miRNA onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Top_miRNA=='') echo 'selected'; echo '>Select Top miRNA</option>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep '.$GMPM.'_ | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "" << "\n";
        output << "        $Count_inFile = File_get_contents( $TSV_List[0] );" << "\n";
        output << "        $Count_inFile_Lines = Explode( \"\\n\", $Count_inFile );" << "\n";
        output << "        $Count_miRNA = Count( $Count_inFile_Lines )-2;" << "\n";
        output << "        $Count_miRNA10 = $Count_miRNA % 10;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 10; $i < $Count_miRNA; $i += 10 )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$i.' ';" << "\n";
        output << "" << "\n";
        output << "            if( $i == $Top_miRNA )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$i.'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<option value='.$Count_miRNA.' ';" << "\n";
        output << "" << "\n";
        output << "        if( $Count_miRNA == $Top_miRNA )" << "\n";
        output << "            echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "        echo '>'.$Count_miRNA.'</option>';" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== mirDistType ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=mirDistType onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $mirDistType_List = array('ppm', '100%');" << "\n";
        output << "        $mirDistType_Size = Count( $mirDistType_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $mirDistType_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$mirDistType_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $mirDistType == $mirDistType_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $mirDistType_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Max NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $mirDistType == 'ppm' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=text name=MaxHight size=8 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $MaxHight=='' )" << "\n";
        output << "                echo 'MaxHight';" << "\n";
        output << "            else" << "\n";
        output << "                echo $MaxHight;" << "\n";
        output << "" << "\n";
        output << "            echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Multi TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File[] multiple=mutiple hight=>';" << "\n";
        output << "" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            foreach( $TSV_File as $TSVs )" << "\n";
        output << "                if( $TSVs == $TSV_List[$i] ) " << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== MirDist ====================-->" << "\n";
        output << "" << "\n";
        output << "        //========== svg view var set ==========" << "\n";
        output << "" << "\n";
        output << "        $width = $Top_miRNA * 180;" << "\n";
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
        output << "        //==================== svg ====================" << "\n";
        output << "" << "\n";
        output << "        if( $mirDistType == 'ppm' )" << "\n";
        output << "        {" << "\n";
        output << "            if( $MaxHight == '' || $MaxHight == 'MaxHight' )" << "\n";
        output << "            {" << "\n";
        output << "                $MaxValue = 0;" << "\n";
        output << "" << "\n";
        output << "                For( $i = 0; $i < Count( $TSV_File ); ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    $inFile = File_get_contents( $TSV_File[$i] );" << "\n";
        output << "                    $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "                    For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "                    {" << "\n";
        output << "                        $inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( $inFile_Line[$k] >= $MaxValue )" << "\n";
        output << "                                $MaxValue = $inFile_Line[$k];" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "                $MaxValue = $MaxHight;" << "\n";
        output << "        }" << "\n";
        output << "        else" << "\n";
        output << "            $MaxValue = 100;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_File ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = File_get_contents( $TSV_File[$i] );" << "\n";
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "            $Temp = Tempnam( '/tmp', $TSV_File[$i] );" << "\n";
        output << "            $Ftemp = fopen( $Temp, 'w' );" << "\n";
        output << "" << "\n";
        output << "            if( $mirDistType == 'ppm' )" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $j >= $Top_miRNA+1 )" << "\n";
        output << "                        break;" << "\n";
        output << "" << "\n";
        output << "                    fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $j >= $Top_miRNA+1 )" << "\n";
        output << "                        break;" << "\n";
        output << "" << "\n";
        output << "                    if( $j == 0 )" << "\n";
        output << "                        fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $inFile_Line_clm = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "                        $Total_Value = 0;" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )" << "\n";
        output << "                            $Total_Value += $inFile_Line_clm[ $k ];" << "\n";
        output << "" << "\n";
        output << "                        fwrite( $Ftemp, $inFile_Line_clm[0] );" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )" << "\n";
        output << "                            fwrite( $Ftemp, \"\\t\".Number_format( $inFile_Line_clm[ $k ]*100/$Total_Value, 0 ));" << "\n";
        output << "" << "\n";
        output << "                        fwrite( $Ftemp, \"\\n\");" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "            echo \"var svg$i = d3.select('body').append('svg')" << "\n";
        output << "                    .attr('id', 'svg$i')" << "\n";
        output << "                .attr('width', width + margin.left + margin.right)" << "\n";
        output << "                .attr('height', height + margin.top + margin.bottom)" << "\n";
        output << "                .append('g')" << "\n";
        output << "                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n";
        output << "            " << "\n";
        output << "                d3.tsv('$Temp', function(error, data) {" << "\n";
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
        output << "                    svg$i.append('g')" << "\n";
        output << "                    .attr('class', 'x axis')" << "\n";
        output << "                    .attr('transform', 'translate(0,' + height + ')')" << "\n";
        output << "                    .call(xAxis);" << "\n";
        output << "            " << "\n";
        output << "                    svg$i.append('g')" << "\n";
        output << "                    .attr('class', 'y axis')" << "\n";
        output << "                    .call(yAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('transform', 'rotate(-90)')" << "\n";
        output << "                    .attr('y', 6)" << "\n";
        output << "                    .attr('dy', '.71em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                        .text('$TSV_File[$i]');" << "\n";
        output << "            " << "\n";
        output << "                    var barchat = svg$i.selectAll('.barchat')" << "\n";
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
        output << "                    var legend = svg$i.selectAll('.legend')" << "\n";
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
