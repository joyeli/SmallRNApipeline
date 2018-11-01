#pragma once
#include <cstdlib>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBoxPlot
{

  public:

    GeneTypeAnalyzerBoxPlot()
    {}

    static void make_file_link( const std::string& p )
    {
        std::string filename;
        std::vector< std::string > list;
        boost::filesystem::path path( p + "../TailDot/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            filename = file.path().filename().string();
            if( filename.substr( filename.length() -3  ) == "tsv" &&
              ( filename.length() < 11 ||
              ( filename.length() > 11 && filename.substr( filename.length() -11 ) != "isomiRs.tsv" ) ))
                list.emplace_back( filename );
        }

        for( auto& file : list )
            if( !boost::filesystem::exists( p + file ))
                 boost::filesystem::create_symlink(( "../TailDot/" + file ).c_str(), ( p + file ).c_str() );
    }

    static void output_boxplot_visualization( const std::string& output_name )
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
        output << "        $ForceY = $_POST['ForceY'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $isTrimmed = $_POST['isTrimmed'];" << "\n";
        output << "        $Filter_Type = $_POST['Filter_Type'];" << "\n";
        output << "        $Filter_Samp = $_POST['Filter_Samp'];" << "\n";
        output << "        $Filter_Pref = $_POST['Filter_Pref'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $TSV_List[$i] != '' )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "                if( $TSV_File == $TSV_List[$i] ) echo 'selected ';" << "\n";
        output << "                echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='hidden' name='Filter_Type' value='$Filter_Type' />" << "\n";
        output << "            <input type='hidden' name='Filter_Samp' value='$Filter_Samp' />" << "\n";
        output << "            <input type='hidden' name='Filter_Pref' value='$Filter_Pref' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================= isTrimmed ===================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isTrimmed onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($isTrimmed=='') echo 'selected'; echo 'value= >isTrimmed</option>';" << "\n";
        output << "" << "\n";
        output << "        $Trim_List = array( '1', '5', '10', '20', '25', '30' );" << "\n";
        output << "        $Trim_Size = Count( $Trim_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Trim_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Trim_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isTrimmed == $Trim_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>±'.$Trim_List[$i].'％</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Filter_Type' value='$Filter_Type' />" << "\n";
        output << "            <input type='hidden' name='Filter_Samp' value='$Filter_Samp' />" << "\n";
        output << "            <input type='hidden' name='Filter_Pref' value='$Filter_Pref' />" << "\n";
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
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='hidden' name='Filter_Type' value='$Filter_Type' />" << "\n";
        output << "            <input type='hidden' name='Filter_Samp' value='$Filter_Samp' />" << "\n";
        output << "            <input type='hidden' name='Filter_Pref' value='$Filter_Pref' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls ../SqAlign/*.idx' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $Sample_List = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_List ) -1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Sample = Explode( '/', $TSV_List[$i] );" << "\n";
        output << "            $Sample = Explode( '.', $Sample[ Count( $Sample ) -1 ] );" << "\n";
        output << "            $Sample_List[ $i ] = $Sample[0];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Type_List = array( 'GMPM', 'GM', 'PM', 'Atail', 'Ctail', 'Gtail', 'Utail' );" << "\n";
        output << "        $Pref_List = array( 'UniLike', 'DisLike' );" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Filter_Samp>';" << "\n";
        output << "        echo '<option value=\"\"'; if($Filter_Samp=='') echo 'selected'; echo '>Filter Sample</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Sample_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Sample_List[$i].' ';" << "\n";
        output << "            if( $Filter_Samp == $Sample_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Sample_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '</select>';" << "\n";
        output << "        echo '<select name=Filter_Type>';" << "\n";
        output << "        echo '<option value=\"\"'; if($Filter_Type=='') echo 'selected'; echo '>Filter Types</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Type_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Type_List[$i].' ';" << "\n";
        output << "            if( $Filter_Type == $Type_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Type_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '</select>';" << "\n";
        output << "        echo '<select name=Filter_Pref>';" << "\n";
        output << "        echo '<option value=\"\"'; if($Filter_Pref=='') echo 'selected'; echo '>Filter Preference</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Pref_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Pref_List[$i].' ';" << "\n";
        output << "            if( $Filter_Pref == $Pref_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Pref_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Read Filter ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Filter_Anno = Array();" << "\n";
        output << "        $Filter_TSV = '../Preference/'.$Filter_Type.'/'.$Filter_Samp.'_'.$Filter_Pref;" << "\n";
        output << "        $Filter_TSV = $Filter_TSV.( SubStr( $TSV_File, 0, 1 ) == 'H' ? '' : '-isomiRs' ).'.text';" << "\n";
        output << "" << "\n";
        output << "        if( $Filter_Samp != '' && $Filter_Type != '' && $Filter_Pref != '' )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = new SplFileObject( $Filter_TSV );" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                if( $inFile_Lines == '' ) continue;" << "\n";
        output << "                Array_Push( $Filter_Anno, Rtrim( $inFile_Lines ));" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "" << "\n";
        output << "        $yMax = 0;" << "\n";
        output << "        $Header = Array();" << "\n";
        output << "" << "\n";
        output << "        $Hete_Array = Array();" << "\n";
        output << "        $Boxs_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File );" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "            $inFile_Lines = Rtrim( $inFile_Lines );" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines );" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                For( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    $Header[ $i-1 ] = $inFile_Line[$i];" << "\n";
        output << "                    $Hete_Array[ $i-1 ] = Array();" << "\n";
        output << "                    $Boxs_Array[ $i-1 ] = Array();" << "\n";
        output << "                }" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            For( $i = 1; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                if( $inFile_Line[$i] == 0 ) continue;" << "\n";
        output << "                if( !Empty( $Filter_Anno ) && !In_Array( $inFile_Line[0], $Filter_Anno )) continue;" << "\n";
        output << "                Array_Push( $Hete_Array[ $i-1 ], $inFile_Line[$i] );" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Hete_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            Sort( $Hete_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "            $Boxs_Array[$i][ 'Q1' ] = $Hete_Array[$i][ Floor(( Count( $Hete_Array[$i] ) -1 ) * 0.25 )];" << "\n";
        output << "            $Boxs_Array[$i][ 'Q2' ] = $Hete_Array[$i][ Floor(( Count( $Hete_Array[$i] ) -1 ) * 0.5  )];" << "\n";
        output << "            $Boxs_Array[$i][ 'Q3' ] = $Hete_Array[$i][ Floor(( Count( $Hete_Array[$i] ) -1 ) * 0.75 )];" << "\n";
        output << "            $Boxs_Array[$i][ 'whisker_low'  ] = $Hete_Array[$i][0];" << "\n";
        output << "            $Boxs_Array[$i][ 'whisker_high' ] = $Hete_Array[$i][ Count( $Hete_Array[$i] ) -1];" << "\n";
        output << "" << "\n";
        output << "            if( $isTrimmed != '' )" << "\n";
        output << "            {" << "\n";
        output << "                $Rate = $isTrimmed / 100;" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_low'  ] = $Hete_Array[$i][ Floor(( Count( $Hete_Array[$i] ) -1 ) *       $Rate  )];" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_high' ] = $Hete_Array[$i][ Floor(( Count( $Hete_Array[$i] ) -1 ) * ( 1 - $Rate ))];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( $Boxs_Array[$i][ 'whisker_high' ] > $yMax ) $yMax = $Boxs_Array[$i][ 'whisker_high' ];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        if( $ForceY != '' ) $yMax = $ForceY;" << "\n";
        output << "" << "\n";
        output << "#<!--================== BoxPlot ====================-->" << "\n";
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
        output << "                var chart = nv.models.boxPlotChart()" << "\n";
        output << "                    .x(function(d) { return d.label })" << "\n";
        output << "                    .staggerLabels(true)" << "\n";
        output << "                    .maxBoxWidth(75) // prevent boxes from being incredibly wide" << "\n";
        output << "                    .yDomain([0, $yMax]);" << "\n";
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
        output << "        For( $i = 0; $i < Count( $Boxs_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '{label:\"'.$Header[$i].'\",values:{';" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Boxs_Array[$i] as $Q => $Value )" << "\n";
        output << "                echo $Q.':'.$Value.',';" << "\n";
        output << "" << "\n";
        output << "            echo 'outliers:[]}}';" << "\n";
        output << "            if( $i != Count( $Boxs_Array )-1 ) echo \",\\n\";" << "\n";
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
