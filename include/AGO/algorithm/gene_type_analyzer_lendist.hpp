#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerLendist
{

  public:

    GeneTypeAnalyzerLendist()
    {}

    static void output_lendist(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::ofstream output( output_name + sample_name + ".tsv" );
        output << "Annotation\tA_Tail\tC_Tail\tG_Tail\tT_Tail\tOther_Tail\tGM";

        for( auto& anno : ano_len_idx.first )
        {
            for( auto& len : ano_len_idx.second )
            {
                output << "\n" << anno << ( anno_mark.find( anno ) != anno_mark.end() ? anno_mark[ anno ] : "" );
                output << ":" << len;

                for( auto& anno_table : anno_table_tail )
                {
                    if( anno_table.find( anno ) != anno_table.end() )
                    {
                        if( anno_table[ anno ].find( len ) != anno_table[ anno ].end() )
                            output << "\t" << std::setprecision( 0 ) << std::fixed << anno_table[ anno ][ len ];
                        else output << "\t0";
                    }
                    else output << "\t0";
                }
            }
        }

        output << "\n";
        output.close();
    }

    static void output_lendist_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "" << "\n";
        output << "        $ForceY = $_POST['ForceY'];" << "\n";
        output << "        $PPMFilter = $_POST['PPMFilter'];" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = $_GET['TSV_File'];" << "\n";
        output << "        $annotation_select = $_GET['annotation_select'];" << "\n";
        output << "" << "\n";
        output << "        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        if( $_GET['annotation_select'] == '' ) $annotation_select = $_POST['annotation_select'];" << "\n";
        output << "        if( $_GET['annotation_select'] == '' ) $isAbundant = $_POST['isAbundant'];" << "\n";
        output << "        else $isAbundant = 'AllAnnotations';" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
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
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== is_Abundant ====================-->" << "\n";
        output << "                " << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n";
        output << "        " << "\n";
        output << "        $isAbundant_List = array('MostAbundant', 'AllAnnotations');" << "\n";
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
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== PPMFilter ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=PPMFilter size=6 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $PPMFilter=='' )" << "\n";
        output << "            echo '\"PPM Filter\"';" << "\n";
        output << "        else" << "\n";
        output << "            echo $PPMFilter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onchange=\\\"form.submit()\\\">\";" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( $TSV_File );" << "\n";
        output << "        $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "        $Annotation_Array = Array();" << "\n";
        output << "        $Anno_Filter = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "            $anno = Explode( ':', $inFile_Line[0] );" << "\n";
        output << "            $Annotation = Explode( '*', $anno[0] );" << "\n";
        output << "" << "\n";
        output << "            if( $isAbundant == 'MostAbundant' && Count( $Annotation ) != 2 )" << "\n";
        output << "                Continue;" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $Anno_Array, $anno[0] );" << "\n";
        output << "" << "\n";
        output << "            if( !isSet( $Anno_Filter[ $anno[0] ] ))" << "\n";
        output << "                $Anno_Filter[ $anno[0] ] = 0;" << "\n";
        output << "" << "\n";
        output << "            for( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                $Anno_Filter[ $anno[0] ] += $inFile_Line[$k];" << "\n";
        output << "" << "\n";
        output << "            if( $annotation_select == $anno[0] || $annotation_select == $Annotation[0] )" << "\n";
        output << "            {" << "\n";
        output << "                $Tail_Array = Array();" << "\n";
        output << "" << "\n";
        output << "                for( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                {" << "\n";
        output << "                    Array_Push( $Tail_Array, $inFile_Line[$k] );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $Annotation_Array[ $anno[1] ] = $Tail_Array;" << "\n";
        output << "                $annotation_select = $anno[0];" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Sort( $Anno_Array, SORT_STRING );" << "\n";
        output << "        $Uniq_Anno_Array = Array_Unique( $Anno_Array, SORT_STRING );" << "\n";
        output << "        $is_any_selected = false;" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=annotation_select onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($annotation_select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Anno_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Uniq_Anno_Array[$i] != '' )" << "\n";
        output << "            {" << "\n";
        output << "                if( $PPMFilter != 'PPM Filter' && $PPMFilter != '' &&" << "\n";
        output << "                    $PPMFilter > $Anno_Filter[ $Uniq_Anno_Array[$i] ])" << "\n";
        output << "                    Continue;" << "\n";
        output << "" << "\n";
        output << "                echo '<option value='.$Uniq_Anno_Array[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $annotation_select == $Uniq_Anno_Array[$i] ) " << "\n";
        output << "                {" << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "                    $is_any_selected = true;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                echo '>'.$Uniq_Anno_Array[$i].'</option>';" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        if( !$is_any_selected ) $Annotation_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== ForceY ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=ForceY size=3 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $ForceY=='' )" << "\n";
        output << "            echo 'Hight';" << "\n";
        output << "        else" << "\n";
        output << "            echo $ForceY;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== GoBack ====================-->" << "\n";
        output << "" << "\n";
        output << "        $annos = Explode( '-', $annotation_select );" << "\n";
        output << "        $anno  = $annos[0];" << "\n";
        output << "" << "\n";
        output << "        For( $i = 1; $i < Count( $annos ) -1; $i++ ) $anno = $anno.'-'.$annos[ $i ];" << "\n";
        output << "" << "\n";
        output << "        if( $TSV_File != '' && $annotation_select != '' )" << "\n";
        output << "            echo \"<a href='../SqAlign/index.php?TSV_File=$TSV_File&Annotation_Select=$anno' >" << "\n";
        output << "                <input type='submit' value='Goto $anno' />" << "\n";
        output << "                </a>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotations Tail Bar Chart ====================-->" << "\n";
        output << "" << "\n";
        output << "        $index = Explode( \"\\t\", $inFile_Lines[0] );" << "\n";
        output << "        $Annotation_Keys = Array_keys( $Annotation_Array );" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $annotation_select );" << "\n";
        output << "        $Ftemp = fopen( $Temp, 'w' );" << "\n";
        output << "        Fwrite( $Ftemp, $index[0] );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Annotation_Keys ); $i++ )" << "\n";
        output << "        {" << "\n";
        output << "            Fwrite( $Ftemp, \"\\t\".$Annotation_Keys[$i] );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $Ftemp, \"\\n\".$index[6] );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Annotation_Keys ); $i++ )" << "\n";
        output << "        {" << "\n";
        output << "            Fwrite( $Ftemp, \"\\t\".$Annotation_Array[ $Annotation_Keys[$i] ][5] );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fwrite( $Ftemp, \"\\n\" );" << "\n";
        output << "" << "\n";
        output << "        For( $j = 1; $j < Count( $index )-1; $j++ )" << "\n";
        output << "        {" << "\n";
        output << "            Fwrite( $Ftemp, $index[$j] );" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $Annotation_Keys ); $i++ )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $Ftemp, \"\\t\".$Annotation_Array[ $Annotation_Keys[$i] ][$j-1] );" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp, \"\\n\" );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "        echo '<svg id=bar></svg>';" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "                d3.tsv( '$Temp', function( tsv_data ) {" << "\n";
        output << "                var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'Annotation'; });" << "\n";
        output << "                var data = new Array();" << "\n";
        output << "        " << "\n";
        output << "                tsv_data.forEach( function(d) {" << "\n";
        output << "                    d.len = length.map( function(key) {" << "\n";
        output << "                        return {" << "\n";
        output << "                            x: Number(key)," << "\n";
        output << "                            y: +d[key]" << "\n";
        output << "                        };" << "\n";
        output << "        " << "\n";
        output << "                    });" << "\n";
        output << "                    var sample = new Object();" << "\n";
        output << "                    sample.key = d['Annotation'];" << "\n";
        output << "                    sample.values = d.len;" << "\n";
        output << "                    data.push( sample );" << "\n";
        output << "                });" << "\n";
        output << "        " << "\n";
        output << "                nv.addGraph({" << "\n";
        output << "                    generate: function() {" << "\n";
        output << "                        var width = nv.utils.windowSize().width," << "\n";
        output << "                            height = nv.utils.windowSize().height;" << "\n";
        output << "            " << "\n";
        output << "                        var chart = nv.models.multiBarChart()\";" << "\n";
        output << "" << "\n";
        output << "            if( $ForceY != '' && $ForceY != 'Hight' )" << "\n";
        output << "                echo \"        .forceY([$ForceY,0])\";" << "\n";
        output << "" << "\n";
        output << "        echo \"                .width(width)" << "\n";
        output << "                            .height(height)" << "\n";
        output << "                            .color( ['#000000', '#FF0000', '#0000FF', '#FFBF00', '#088A08', '#6E6E6E'])" << "\n";
        output << "                            .stacked(true);" << "\n";
        output << "        " << "\n";
        output << "                        var svg = d3.select('#bar').attr('height',height*2).datum(data);" << "\n";
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
        output << "                                d3.select('#bar')" << "\n";
        output << "                                        .attr('width', width)" << "\n";
        output << "                                        .attr('height', height)" << "\n";
        output << "                                        .transition().duration(0)" << "\n";
        output << "                                        .call(graph);" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "                });" << "\n";
        output << "            });" << "\n";
        output << "        </script>\";" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
