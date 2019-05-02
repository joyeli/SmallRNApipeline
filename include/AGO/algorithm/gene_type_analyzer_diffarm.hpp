#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDiffarm
{
  public:

    using DiffArmType = std::map< std::string, std::vector< std::vector< double >>>;
    //                              anno       5p[0] 3p[1]  Tail%[0] A[1] C[2] G[3] U[4] O[5]

    GeneTypeAnalyzerDiffarm()
    {}

    static void make_file_link( const std::string& p )
    {
        std::string filename;
        std::vector< std::string > list;
        std::vector< std::string > split;

        boost::filesystem::path path( p + "../TailDot/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            filename = file.path().filename().string();

            boost::iter_split( split, filename, boost::algorithm::first_finder( "-isomiRs" ));
            if( split.size() > 1 ) continue;

            if( filename.substr( filename.length() -3  ) == "tsv" )
                list.emplace_back( filename );
        }

        for( auto& file : list )
            if( !boost::filesystem::exists( p + file ))
                 boost::filesystem::create_symlink(( "../TailDot/" + file ).c_str(), ( p + file ).c_str() );
    }

    static void output_diffarm_visualization( const std::string& output_name )
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
        output << "        $Type = $_POST['Type'];" << "\n";
        output << "        $isPair = $_POST['isPair'];" << "\n";
        output << "        $Sample = $_POST['Sample'];" << "\n";
        output << "        $FilterTop = $_POST['FilterTop'];" << "\n";
        output << "        $FilterDwn = $_POST['FilterDwn'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== GetSamp =====================-->" << "\n";
        output << "        " << "\n";
        output << "        $Samp = Shell_Exec( 'ls *.tsv' );" << "\n";
        output << "        $Samp_List = Explode( \"\\n\", $Samp );" << "\n";
        output << "        $Sample_List = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Samp_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Samp_List[$i] == '' ) continue;" << "\n";
        output << "            $Samp_Name = Explode( '.tsv', $Samp_List[$i] );" << "\n";
        output << "            $Samp_Name = Explode( '/', $Samp_Name[ 0] );" << "\n";
        output << "            Array_Push( $Sample_List, $Samp_Name[ Count( $Samp_Name ) -1 ]);" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== isPair =====================-->" << "\n";
        output << "                " << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=isPair onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isPair=='') echo 'selected'; echo '>isPair</option>';" << "\n";
        output << "" << "\n";
        output << "        $List = array('Yes', 'No');" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $List ); ++$i )" << "\n";
        output << "            echo '<option value='.$List[$i].( $isPair == $List[$i] ? ' selected >' : ' >' ).$List[$i].'</option>';" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='Sample' value='$Sample' />" << "\n";
        output << "            <input type='hidden' name='FilterTop' value='$FilterTop' />" << "\n";
        output << "            <input type='hidden' name='FilterDwn' value='$FilterDwn' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Sample =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Sample onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Sample=='') echo 'selected'; echo '>Select Sample</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Sample_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Sample_List[$i] == '' ) continue;" << "\n";
        output << "            echo '<option value='.$Sample_List[$i].' ';" << "\n";
        output << "            if( $Sample == $Sample_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Sample_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='FilterTop' value='$FilterTop' />" << "\n";
        output << "            <input type='hidden' name='FilterDwn' value='$FilterDwn' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=============== Type ==================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Type onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Type=='') echo 'selected'; echo '>Select Types</option>';" << "\n";
        output << "" << "\n";
        output << "        $inFileHeader = new SplFileObject( $Sample.'.tsv' );" << "\n";
        output << "        $inHeaders = Explode( \"\\t\", $inFileHeader->fgets() );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 1; $i < Count( $inHeaders ); ++$i )" << "\n";
        output << "            echo '<option value='.$i.( $Type == $i ? ' selected' : '' ).' >'.$inHeaders[$i].'</option>';" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='Sample' value='$Sample' />" << "\n";
        output << "            <input type='hidden' name='FilterTop' value='$FilterTop' />" << "\n";
        output << "            <input type='hidden' name='FilterDwn' value='$FilterDwn' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== FilterGMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=FilterTop size=3 value=';" << "\n";
        output << "        echo $FilterTop == '' ? 'FilterTop' : $FilterTop;" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=FilterDwn size=3 value=';" << "\n";
        output << "        echo $FilterDwn == '' ? 'FilterDwn' : $FilterDwn;" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='Type' value='$Type' />" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='Sample' value='$Sample' />" << "\n";
        output << "            <input type='submit' value='Submit' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "" << "\n";
        output << "        $Value_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $inFile = new SplFileObject( $Sample.'.tsv' );" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "            if( $inFile_Line[0] == $Sample ) continue;" << "\n";
        output << "            if( $inFile_Line[1] == 0 ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Line[0] = Explode( '!', $inFile_Line[0] )[0];" << "\n";
        output << "            $miRNA_Arm = Explode( '-', $inFile_Line[0] );" << "\n";
        output << "" << "\n";
        output << "            $Arm = $miRNA_Arm[ Count( $miRNA_Arm ) -1 ];" << "\n";
        output << "            $miRNA = $miRNA_Arm[0];" << "\n";
        output << "" << "\n";
        output << "            for( $i = 1 ; $i < Count( $miRNA_Arm ) -1; ++$i )" << "\n";
        output << "                $miRNA = $miRNA.'-'.$miRNA_Arm[ $i ];" << "\n";
        output << "" << "\n";
        output << "            if( !array_key_exists( $miRNA, $Value_Array ))" << "\n";
        output << "            {" << "\n";
        output << "                $Value_Array[ $miRNA ] = Array();" << "\n";
        output << "                $Value_Array[ $miRNA ][ 'GMPM' ] = 0;" << "\n";
        output << "                $Value_Array[ $miRNA ][ 'SUM' ] = 0;" << "\n";
        output << "                $Value_Array[ $miRNA ][ '5p' ] = 0;" << "\n";
        output << "                $Value_Array[ $miRNA ][ '3p' ] = 0;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Value_Array[ $miRNA ][ 'GMPM' ] += $inFile_Line[1];" << "\n";
        output << "            $Value_Array[ $miRNA ][ $Arm ] = $inFile_Line[ $Type ];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== SumValues =====================-->" << "\n";
        output << "" << "\n";
        output << "        foreach( $Value_Array as $miR => $Value )" << "\n";
        output << "        {" << "\n";
        output << "            if(( $Value_Array[ $miR ][ 'GMPM' ] == 0 ) ||" << "\n";
        output << "               ( $FilterTop != '' && $FilterTop != 'FilterTop' && Abs( $Value_Array[ $miR ][ 'GMPM' ]) > $FilterTop ) ||" << "\n";
        output << "               ( $FilterDwn != '' && $FilterDwn != 'FilterDwn' && Abs( $Value_Array[ $miR ][ 'GMPM' ]) < $FilterDwn ) ||" << "\n";
        output << "               ( $isPair == 'Yes' && ( $Value_Array[ $miR ][ '5p' ] == 0 || $Value_Array[ $miR ][ '3p' ] == 0 )     ) ||" << "\n";
        output << "               ( $Value_Array[ $miR ][ '5p' ] == 0 && $Value_Array[ $miR ][ '3p' ] == 0 ))" << "\n";
        output << "            {" << "\n";
        output << "                unset( $Value_Array[ $miR ] );" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Value_Array[ $miR ][ 'SUM' ] = $Value_Array[ $miR ][ '5p' ] + $Value_Array[ $miR ][ '3p' ];" << "\n";
        output << "" << "\n";
        output << "            if( $Value_Array[ $miR ][ 'SUM' ] != 0 )" << "\n";
        output << "            {" << "\n";
        output << "                $Value_Array[ $miR ][ '5p' ] = $Value_Array[ $miR ][ '5p' ] / $Value_Array[ $miR ][ 'SUM' ];" << "\n";
        output << "                $Value_Array[ $miR ][ '3p' ] = $Value_Array[ $miR ][ '3p' ] / $Value_Array[ $miR ][ 'SUM' ];" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Sort_Array = Array();" << "\n";
        output << "        $Data_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $Data_Array[ '5p' ] = Array();" << "\n";
        output << "        $Data_Array[ '3p' ] = Array();" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Value_Array as $miR => $Value )" << "\n";
        output << "            Array_Push( $Sort_Array, $Value_Array[ $miR ][ '5p' ] - $Value_Array[ $miR ][ '3p' ]);" << "\n";
        output << "" << "\n";
        output << "        Array_Multisort( $Sort_Array, SORT_NUMERIC, SORT_DESC, $Value_Array );" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Value_Array as $miR => $Value )" << "\n";
        output << "        {" << "\n";
        output << "            $Data_Array[ '5p' ][ $miR.' : '.$Value_Array[ $miR ][ 'GMPM' ]] = $Value_Array[ $miR ][ '5p' ];" << "\n";
        output << "            $Data_Array[ '3p' ][ $miR.' : '.$Value_Array[ $miR ][ 'GMPM' ]] = $Value_Array[ $miR ][ '3p' ] * -1;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== HisGram ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "            var margin = {top: 10, right: 20, bottom: 10, left: 20}" << "\n";
        output << "                width  = svg_width  - margin.left - margin.right," << "\n";
        output << "                height = svg_height - margin.top  - margin.bottom;" << "\n";
        output << "" << "\n";
        output << "            var svg = d3.select('body').append('div')" << "\n";
        output << "                .attr('id', 'svg')" << "\n";
        output << "                .style('width', width + 'px')" << "\n";
        output << "                .style('height', height + 'px')" << "\n";
        output << "                .style('display', 'inline-block' )" << "\n";
        output << "                .append('svg');" << "\n";
        output << "" << "\n";
        output << "            nv.addGraph( function() {" << "\n";
        output << "                var chart = nv.models.multiBarHorizontalChart()" << "\n";
        output << "                    .x(function(d) { return d.label })" << "\n";
        output << "                    .y(function(d) { return d.value })" << "\n";
        output << "                    .color([ '#d67777', '#4f99b4' ])" << "\n";
        output << "                    .margin({left: 140})" << "\n";
        output << "                    .groupSpacing(0)" << "\n";
        output << "                    .showControls(true)" << "\n";
        output << "                    .showValues(true)" << "\n";
        output << "                    .stacked(true)" << "\n";
        output << "                    .forceY([-1,0,1])" << "\n";
        output << "                    ;" << "\n";
        output << "" << "\n";
        output << "                d3.select('#svg svg')" << "\n";
        output << "                    .data([data])" << "\n";
        output << "                    .call(chart);" << "\n";
        output << "" << "\n";
        output << "                nv.utils.windowResize(chart.update);" << "\n";
        output << "                return chart;" << "\n";
        output << "            });" << "\n";
        output << "" << "\n";
        output << "            var data = [\";" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Data_Array as $Key => $Data )" << "\n";
        output << "        {" << "\n";
        output << "            echo '{key:\"'.$Key.'\",values:[';" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Data as $Label => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                echo '{label:\"'.$Label.'\",value:'.$Value.'},';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo ']},';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '];</script>';" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
