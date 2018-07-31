#pragma once

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerSeedPie
{
  public:

    GeneTypeAnalyzerSeedPie()
    {}

    static void make_file_link( const std::string& p )
    {
        std::vector< std::string > list;
        boost::filesystem::path path( p + "../SeedMap/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            std::string name = file.path().filename().string();
            if( name.substr( name.length() -3, 3 ) == "php" ) continue;
            list.emplace_back( name );
        }

        for( auto& file : list )
            if( !boost::filesystem::exists( p + file ))
                 boost::filesystem::create_symlink(( "../SeedMap/" + file ).c_str(), ( p + file ).c_str() );
    }

    static void output_seedpie_visualization( const std::string& output_name )
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
        output << "" << "\n";
        output << "        $Filter = $_GET['Filter'];" << "\n";
        output << "        $TSV_File = $_GET['TSV_File'];" << "\n";
        output << "        $Annotation_Select = $_GET['Annotation_Select'];" << "\n";
        output << "" << "\n";
        output << "        if( $_GET['Filter'] == '' ) $Filter = $_POST['Filter'];" << "\n";
        output << "        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        if( $_GET['Annotation_Select'] == '' ) $Annotation_Select = $_POST['Annotation_Select'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js ></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
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
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "        $Total_PPM = 0.0;" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $PPM = 0.0;" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            For( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                $PPM += $inFile_Line[$i];" << "\n";
        output << "" << "\n";
        output << "            $Anno_Array[ $inFile_Line[0] ] = $PPM;" << "\n";
        output << "            if( $PPM > $Total_PPM ) $Total_PPM = $PPM;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Annotation_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$anno.' ';" << "\n";
        output << "            if( $Annotation_Select == $anno ) echo 'selected ';" << "\n";
        output << "            echo '>'.$anno.' ('.$ppm.'ppm)</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "        $Total_PPM = 100 / $Total_PPM;" << "\n";
        output << "        echo \"<script>var select_color_map = d3.scale.linear().domain([ 0, 100 ]).range([ 'WhiteSmoke', 'Black' ]);\";" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "        {" << "\n";
        output << "           echo \"$( 'option[value=\\\"$anno\\\"]' ).css(" << "\n";
        output << "           {" << "\n";
        output << "               'background-color': select_color_map( '\".( $ppm * $Total_PPM ).\"' )," << "\n";
        output << "               'color': '\".(( $ppm * $Total_PPM ) > 25 ? 'white' : 'black'  ).\"'" << "\n";
        output << "           });\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '</script>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == '' ) $Filter = 'FilterGMPM';" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value='.$Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== GoLenD ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $TSV_File != '' && $Annotation_Select != '' && $Annotation_Select != 'Select Annotations' )" << "\n";
        output << "        {" << "\n";
        output << "            $TSV = Substr( $TSV_File, 0, Strlen( $TSV_File ) -4 );" << "\n";
        output << "            echo \"<a target='_blank' href='../LenDist/index.php?TSV_File=$TSV&annotation_select=$Annotation_Select' >" << "\n";
        output << "                <input type='submit' value='Goto $Annotation_Select' />" << "\n";
        output << "                </a>\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--===================== Read File ======================-->" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $TSV_File.'_'.$Annotation_Select );" << "\n";
        output << "        $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                for( $i = 1; $i < Count( $inFile_Line ); ++$i ) $Anno_Array[$i] = Rtrim( $inFile_Line[$i] );" << "\n";
        output << "                Fwrite( $Ftemp, \"Anno\\tPPM\\n\" );" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( $inFile_Line[0] != $Annotation_Select ) continue;" << "\n";
        output << "" << "\n";
        output << "            for( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Line[$i] = Rtrim( $inFile_Line[$i]);" << "\n";
        output << "                if( $Filter != '' && $Filter != 'FilterGMPM' && $inFile_Line[$i] < $Filter ) continue;" << "\n";
        output << "                Fwrite( $Ftemp, $Anno_Array[$i].\"\\t\".$inFile_Line[$i].\"\\n\" );" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            break;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "#<!--================== PieChart ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<svg id=pie></svg>';" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "                d3.tsv( '$Temp', function( tsv_data ) {" << "\n";
        output << "" << "\n";
        output << "                var data = new Array();" << "\n";
        output << "        " << "\n";
        output << "                tsv_data.forEach( function(d) {" << "\n";
        output << "                    var annotation = new Object();" << "\n";
        output << "                    annotation.key = d['Anno'];" << "\n";
        output << "                    annotation.value = d['PPM'];" << "\n";
        output << "                    data.push( annotation );" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "                nv.addGraph({" << "\n";
        output << "                    generate: function() {" << "\n";
        output << "                        var width = nv.utils.windowSize().width," << "\n";
        output << "                            height = nv.utils.windowSize().height;" << "\n";
        output << "            " << "\n";
        output << "                        var chart = nv.models.pieChart()" << "\n";
        output << "                            .width(width)" << "\n";
        output << "                            .height(height)" << "\n";
        output << "                            .x(function(d) { return d.key })" << "\n";
        output << "                            .y(function(d) { return d.value })" << "\n";
        output << "                            .showLabels(true)     //Display pie labels" << "\n";
        output << "                            .labelThreshold(0.009)//Configure the minimum slice size for labels to show up" << "\n";
        output << "                            .labelType('key') //Configure what type of data to show in the label. Can be 'key', 'value' or 'percent'" << "\n";
        output << "                            .labelsOutside(true)" << "\n";
        output << "        " << "\n";
        output << "                        var svg = d3.select('#pie').attr('height',height).datum(data);" << "\n";
        output << "                        svg.transition().duration(0).call(chart);" << "\n";
        output << "            " << "\n";
        output << "                        return chart;" << "\n";
        output << "                    }," << "\n";
        output << "                    callback: function(graph) {" << "\n";
        output << "                            nv.utils.windowResize(function() {" << "\n";
        output << "                                var width = nv.utils.windowSize().width;" << "\n";
        output << "                                var height = nv.utils.windowSize().height;" << "\n";
        output << "                                graph.width(width).height(height);" << "\n";
        output << "" << "\n";
        output << "                                d3.select('#pie')" << "\n";
        output << "                                        .attr('width', width)" << "\n";
        output << "                                        .attr('height', height)" << "\n";
        output << "                                        .transition().duration(0)" << "\n";
        output << "                                        .call(graph);" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "                });" << "\n";
        output << "            });" << "\n";
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
