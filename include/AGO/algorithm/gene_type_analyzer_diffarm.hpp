#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerDiffarm
{
  public:

    using DiffArmType = std::map< std::string, std::vector< std::vector< double >>>;
    //                              anno       5p[0] 3p[1]  Tail%[0] A[1] C[2] G[3] T[4] O[5]

    GeneTypeAnalyzerDiffarm()
    {}

    static std::pair< std::string, std::string > get_anno_arm( const std::string& anno )
    {
        std::vector< std::string > split;
        std::pair< std::string, std::string > anno_arm;

        boost::iter_split( split, anno, boost::algorithm::first_finder( "-" ));
        std::string name = split[0];

        for( std::size_t i = 1; i < split.size() -1; ++i )
            name += "-" + split[i];

        return { name, split[ split.size() -1 ].substr( 0, 2 ) };
    }

    static DiffArmType read_tsv( const std::string& file )
    {
        DiffArmType diffarm;
        std::string line;

        std::ifstream infile( file );
        std::vector< std::string > split;

        std::pair< std::string, std::string > anno_arm;

        while( std::getline( infile, line ))
        {
            boost::iter_split( split, line, boost::algorithm::first_finder( "\t" ));
            if( split[1] == "GMPM" ) continue;

            anno_arm = get_anno_arm( split[0] );
            if( diffarm.find( anno_arm.first ) == diffarm.end() )
                diffarm[ anno_arm.first ] = std::vector< std::vector< double >>( 2, std::vector< double >( 8, 0.000001 ));

            diffarm[ anno_arm.first ][ anno_arm.second == "5p" ? 0 : 1 ] = 
            {
                std::stod( split[1] ),
                std::stod( split[1] ) * std::stod( split[3] ),
                std::stod( split[3] ) + 0.000001,
                std::stod( split[4] ) + 0.000001,
                std::stod( split[5] ) + 0.000001,
                std::stod( split[6] ) + 0.000001,
                std::stod( split[7] ) + 0.000001,
                std::stod( split[8] ) + 0.000001
            };
        }

        infile.close();
        return diffarm;
    }

    static void get_diff( DiffArmType& diffarm )
    {
        for( auto& anno : diffarm )
        {
            auto& p5 = anno.second[0];
            auto& p3 = anno.second[1];

            p5[0] = ( p5[0] == 0.000001 ? 0 : p5[0] ) + ( p3[0] == 0.000001 ? 0 : p3[0] );
            p5[1] =   p5[0] == 0 ? 0 : ((( p5[1] == 0.000001 ? 0 : p5[1] ) + ( p3[1] == 0.000001 ? 0 : p3[1] )) / p5[0] );

            for( std::size_t i = 2; i < 8; ++i )
            {
                auto& p5 = anno.second[0][i];
                auto& p3 = anno.second[1][i];
                p5 = p5 - p3;
                //p5 = p5 == p3 && p5 == 0.000001 ? 0 : (( p5 > p3 ? p5 : p3 ) / (( p5 == 0.000001 ? 0 : p5 ) + ( p3 == 0.000001 ? 0 : p3 )));
            }
        }
    }

    static void out_diff( DiffArmType& diffarm, const std::string& out_file )
    {
        std::ofstream output( out_file );
        output << "Anno\tGMPM\tTailingRatio\tDiff_Tailing\tDiff_A\tDiff_C\tDiff_G\tDiff_T\tDiff_Other\n";

        for( auto& anno : diffarm )
        {
            output << anno.first;

            for( auto& ratio : anno.second[0] )
                output << "\t" << ratio;

            output << "\n";
        }

        output.close();
    }

    static void get_diffarm( const std::string& p )
    {
        std::string filename;
        std::vector< std::string > list;
        std::vector< std::string > split;

        boost::filesystem::path path( p + "../TailDot/" );
        DiffArmType diffarm;

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            filename = file.path().filename().string();

            boost::iter_split( split, filename, boost::algorithm::first_finder( "-isomiRs" ));
            if( split.size() != 1 ) continue;

            boost::iter_split( split, filename, boost::algorithm::first_finder( ".tsv" ));
            if( split.size() == 1 ) continue;

            if( filename.substr( 0, 15 ) != "Heterorgeneity_" )
                list.emplace_back( filename );
        }

        for( auto& file : list )
        {
            diffarm = read_tsv( p + "../TailDot/" + file );
            get_diff( diffarm );
            out_diff( diffarm, p + file );
        }
    }

    static void output_diffarm_visualization( const std::string& output_name )
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
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $DiffType = $_POST['DiffType'];" << "\n";
        output << "        $isTrimmed = $_POST['isTrimmed'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Get Header ====================-->" << "\n";
        output << "        " << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $TSVs = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_List ); ++$i )" << "\n";
        output << "            if( $TSV_List[$i] != '' ) Array_Push( $TSVs, $TSV_List[$i] );" << "\n";
        output << "" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $inFile = new SplFileObject( $TSVs[0] );" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "" << "\n";
        output << "                For( $i = 3; $i < Count( $inFile_Line ); ++$i )" << "\n";
        output << "                    $Header[ $inFile_Line[$i] ] = $i;" << "\n";
        output << "" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "            else break;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=DiffType onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($DiffType=='') echo 'selected'; echo 'value= >DiffType</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Header as $Type => $Index )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Type.' ';" << "\n";
        output << "" << "\n";
        output << "            if( $DiffType == $Type )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Type.'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
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
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == '' ) $Filter = 'FilterGMPM';" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value='.$Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "" << "\n";
        output << "        $yMin = 0;" << "\n";
        output << "        $yMax = 0;" << "\n";
        output << "" << "\n";
        output << "        $Tail_Array = Array();" << "\n";
        output << "        $Boxs_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        if( $DiffType != '' && $DiffType != 'DiffType' )" << "\n";
        output << "            For( $i = 0; $i < Count( $TSVs ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $isHeader = true;" << "\n";
        output << "                $Tail_Array[$i] = Array();" << "\n";
        output << "                $inFile = new SplFileObject( $TSVs[$i] );" << "\n";
        output << "" << "\n";
        output << "                while( !$inFile->eof() )" << "\n";
        output << "                {" << "\n";
        output << "                    $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                    if( $inFile_Lines == '' ) continue;" << "\n";
        output << "                    $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "                    if( $isHeader )" << "\n";
        output << "                    {" << "\n";
        output << "                        $isHeader = false;" << "\n";
        output << "                        continue;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( $inFile_Line[2] == 0 ) continue;" << "\n";
        output << "                    if( $Filter != '' && $Filter != 'FilterGMPM' && $inFile_Line[1] < $Filter ) continue;" << "\n";
        output << "                    Array_Push( $Tail_Array[$i], $inFile_Line[ $Header[ $DiffType ]] );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Tail_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            Sort( $Tail_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "            $Boxs_Array[$i][ 'Q1' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.25 )];" << "\n";
        output << "            $Boxs_Array[$i][ 'Q2' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.5  )];" << "\n";
        output << "            $Boxs_Array[$i][ 'Q3' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.75 )];" << "\n";
        output << "            $Boxs_Array[$i][ 'whisker_low'  ] = $Tail_Array[$i][0];" << "\n";
        output << "            $Boxs_Array[$i][ 'whisker_high' ] = $Tail_Array[$i][ Count( $Tail_Array[$i] ) -1];" << "\n";
        output << "" << "\n";
        output << "            if( $isTrimmed != '' && $isTrimmed != 'isTrimmed' )" << "\n";
        output << "            {" << "\n";
        output << "                $Rate = $isTrimmed / 100;" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_low'  ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) *       $Rate  )];" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_high' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * ( 1 - $Rate ))];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( $Boxs_Array[$i][ 'whisker_low'  ] < $yMin ) $yMin = $Boxs_Array[$i][ 'whisker_low'  ];" << "\n";
        output << "            if( $Boxs_Array[$i][ 'whisker_high' ] > $yMax ) $yMax = $Boxs_Array[$i][ 'whisker_high' ];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== BoxPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        $yMin = -1;" << "\n";
        output << "        $yMax = 1;" << "\n";
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
        output << "                    .yDomain([$yMin, $yMax]);" << "\n";
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
        output << "            echo '{label:\"'.Explode( '.', $TSVs[$i] )[0].'\",values:{';" << "\n";
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
