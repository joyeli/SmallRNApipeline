#pragma once

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerHisGram
{
  public:

    GeneTypeAnalyzerHisGram()
    {}

    static std::map< std::string, std::vector< double >> read_tsv( const std::string& output_name, std::map< std::size_t, std::string >& headers )
    {
        bool is_skip = false;
        std::size_t line_count = 0;
        std::string line;
        std::ifstream infile( output_name );

        std::vector< double > values;
        std::vector< std::string > split;
        std::map< std::string, std::vector< double >> table;

        while( std::getline( infile, line ))
        {
            line_count++;
            values.clear();
            is_skip = false;

            boost::iter_split( split, line, boost::algorithm::first_finder( "\t" ));

            for( std::size_t i = 1; i < split.size(); ++i )
            {
                if( line_count != 1 )
                {
                    if( split[i] == "0" )
                    {
                        is_skip = true;
                        break;
                    }
                    else values.emplace_back( std::stod( split[i] ));
                }
                else
                {
                    is_skip = true;
                    headers.emplace( i-1, split[i] );
                }
            }

            if( is_skip ) continue;
            table[ split[0] ] = values;
        }

        return table;
    }

    static std::vector< std::vector< std::pair< std::size_t, std::string >>>
        get_smp_compare( std::map< std::size_t, std::string >& headers )
    {
        std::vector< std::pair< std::size_t, std::string >> compares_temp;
        std::vector< std::vector< std::pair< std::size_t, std::string >>> compares;

        std::set< std::pair< std::size_t, std::string >> smp_check_temp;
        std::set< std::set< std::pair< std::size_t, std::string >>> smp_check;

        for( auto& header1 : headers )
        {
            for( auto& header2 : headers )
            {
                if( header1.second == header2.second ) continue;
                smp_check_temp.emplace( header1 );
                smp_check_temp.emplace( header2 );
                smp_check.emplace( smp_check_temp );
                smp_check_temp.clear();
            }
        }

        for( auto& compare : smp_check )
        {
            for( auto& smp : compare ) compares_temp.emplace_back( smp );
            compares.emplace_back( compares_temp );
            compares_temp.clear();
        }

        return compares;
    }

    static int get_bin( const int& value, double bin )
    {
        int tmp = (int)( value < 0 ? (-bin) : bin );
        int res = (( value + tmp ) - ( value + tmp ) % tmp );
        return ( bin < 0 ? ( res - tmp ) : res );
    }

    static std::map< double, std::size_t > make_hisgram( auto& dists, double bin = 0.025, int shift = 1000 )
    {
        int min = 0;
        int max = 0;
        int tmp = 0;

        bin = bin * shift;
        std::map< int, std::size_t > hisgrams_temp;
        std::map< double, std::size_t > hisgrams;

        for( auto dist : dists )
        {
            tmp = dist.second * shift;
            if( tmp < min ) min = tmp;
            if( tmp > max ) max = tmp;
        }

        min = get_bin( min, bin );
        max = get_bin( max, bin );

        for( int i = min; i <= max; i += bin )
        {
            if( i > ((int)-bin) && i < (int)bin ) i = 0;
            hisgrams_temp[i] = 0;
        }

        for( auto& dist : dists )
        {
            tmp = dist.second * shift;
            tmp = get_bin( tmp, bin );
            hisgrams_temp[ tmp ]++;
        }

        for( auto& histram : hisgrams_temp )
            hisgrams[ (double)histram.first / (double)shift ] = histram.second;

        return hisgrams;
    }

    static void output_hisgram_from_heter( const std::string& output_name )
    {
        if( !boost::filesystem::exists( output_name + "../BoxPlot/Heterorgeneity_5p.tsv" ) ||
            !boost::filesystem::exists( output_name + "../BoxPlot/Heterorgeneity_3p.tsv" ) )
            return;

        std::map< std::string, double > dists;
        std::map< double, std::size_t > hisgrams;
        std::map< std::size_t, std::string > headers;
        std::map< std::string, std::map< std::string, std::vector< double >>> tables;

        std::vector< std::vector< std::pair< std::size_t, std::string >>> compares;
        std::ofstream output;

        tables[ "5p" ] = read_tsv( output_name + "../BoxPlot/Heterorgeneity_5p.tsv", headers );
        tables[ "3p" ] = read_tsv( output_name + "../BoxPlot/Heterorgeneity_3p.tsv", headers );

        for( auto& compare : get_smp_compare( headers ))
        {
            auto& s1 = compare[0];
            auto& s2 = compare[1];

            for( auto& table : tables )
            {
                for( auto& anno : table.second )
                    dists[ anno.first ] = anno.second[ s1.first ] - anno.second[ s2.first ];

                output.open( output_name + "Heterorgeneity_" + table.first + "_" + s1.second + "_" + s2.second + "_raw.tsv" );
                output << "Anno\t" << s1.second << "-" << s2.second;

                for( auto& dist : dists )
                    output << "\n" << dist.first << "\t" << dist.second;

                output.close();

                hisgrams = make_hisgram( dists );

                output.open( output_name + "Heterorgeneity_" + table.first + "_" + s1.second + "_" + s2.second + ".tsv" );
                output << "Bin\t" << s1.second << "-" << s2.second;

                for( auto& hisgram : hisgrams )
                    output << "\n" << hisgram.first << "\t" << hisgram.second;

                output.close();
                hisgrams.clear();
                dists.clear();
            }
        }
    }

    static void output_hisgram_from_entro( const std::string& output_name )
    {
        if( !boost::filesystem::exists( output_name + "../BoxPlot/Entropy_5p.tsv"  ) ||
            !boost::filesystem::exists( output_name + "../BoxPlot/Entropy_mid.tsv" ) ||
            !boost::filesystem::exists( output_name + "../BoxPlot/Entropy_3p.tsv"  ) ||
            !boost::filesystem::exists( output_name + "../BoxPlot/Entropy_3pTailOnly.tsv" ))
            return;

        std::map< std::string, double > dists;
        std::map< double, std::size_t > hisgrams;
        std::map< std::size_t, std::string > headers;
        std::map< std::string, std::map< std::string, std::vector< double >>> tables;

        std::vector< std::vector< std::pair< std::size_t, std::string >>> compares;
        std::ofstream output;

        tables[ "5p"  ] = read_tsv( output_name + "../BoxPlot/Entropy_5p.tsv" , headers );
        tables[ "mid" ] = read_tsv( output_name + "../BoxPlot/Entropy_mid.tsv", headers );
        tables[ "3p"  ] = read_tsv( output_name + "../BoxPlot/Entropy_3p.tsv" , headers );
        tables[ "3pTailOnly" ] = read_tsv( output_name + "../BoxPlot/Entropy_3pTailOnly.tsv", headers );

        for( auto& compare : get_smp_compare( headers ))
        {
            auto& s1 = compare[0];
            auto& s2 = compare[1];

            for( auto& table : tables )
            {
                for( auto& anno : table.second )
                    dists[ anno.first ] = anno.second[ s1.first ] - anno.second[ s2.first ];

                output.open( output_name + "Entropy_" + table.first + "_" + s1.second + "_" + s2.second + "_raw.tsv" );
                output << "Anno\t" << s1.second << "-" << s2.second;

                for( auto& dist : dists )
                    output << "\n" << dist.first << "\t" << dist.second;

                output.close();

                hisgrams = make_hisgram( dists );

                output.open( output_name + "Entropy_" + table.first + "_" + s1.second + "_" + s2.second + ".tsv" );
                output << "Bin\t" << s1.second << "-" << s2.second;

                for( auto& hisgram : hisgrams )
                    output << "\n" << hisgram.first << "\t" << hisgram.second;

                output.close();
                hisgrams.clear();
                dists.clear();
            }
        }
    }

    static void output_hisgram_visualization( const std::string& output_name )
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
        output << "        $Bin_Size = $_POST['Bin_Size'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $TSV_Sample = $_POST['TSV_Sample'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep _raw.tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $TSV_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $TSV_List[$i] == '' ) continue;" << "\n";
        output << "            $TSV_Name = Explode( '_', $TSV_List[$i] );" << "\n";
        output << "            $TSV_Temp = $TSV_Name[0].'_'.$TSV_Name[1];" << "\n";
        output << "" << "\n";
        output << "            if( !Array_Key_Exists( $TSV_Temp, $TSV_Array ))" << "\n";
        output << "                $TSV_Array[ $TSV_Temp ] = Array();" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $TSV_Array[ $TSV_Temp ], $TSV_Name[2].'_'.$TSV_Name[3] );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo 'value= >Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $TSV_Array as $TSV_Name => $TSV_Files )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_Name.' ';" << "\n";
        output << "            if( $TSV_File == $TSV_Name ) echo 'selected ';" << "\n";
        output << "            echo '>'.$TSV_Name.'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Bin_Size' value='$Bin_Size' />" << "\n";
        output << "            <input type='hidden' name='TSV_Sample' value='$TSV_Sample' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV Sample ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=TSV_Sample onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_Sample=='') echo 'selected'; echo 'value= >Select Sample</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $TSV_Array[ $TSV_File ] as $Samples )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Samples.' ';" << "\n";
        output << "            if( $TSV_Sample == $Samples ) echo 'selected ';" << "\n";
        output << "            echo '>'.Str_Replace( '_', '-', $Samples ).'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Bin_Size' value='$Bin_Size' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Bin =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Bin_Size onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Bin_Size=='') echo 'selected'; echo 'value= >Bin Size</option>';" << "\n";
        output << "" << "\n";
        output << "        $Bin_List = array( 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1 );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Bin_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Bin_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Bin_Size == $Bin_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Bin_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='TSV_Sample' value='$TSV_Sample' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "" << "\n";
        output << "        $Value_Array = Array();" << "\n";
        output << "        $Min = 0.0;" << "\n";
        output << "        $Max = 0.0;" << "\n";
        output << "" << "\n";
        output << "        $Header = '';" << "\n";
        output << "        $isHeader = true;" << "\n";
        output << "        $inFile = new SplFileObject( $TSV_File.'_'.$TSV_Sample.'_raw.tsv' );" << "\n";
        output << "" << "\n";
        output << "        while( !$inFile->eof() )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Lines = $inFile->fgets();" << "\n";
        output << "            if( $inFile_Lines == '' ) continue;" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "" << "\n";
        output << "            if( $isHeader )" << "\n";
        output << "            {" << "\n";
        output << "                $Header = $inFile_Line[1];" << "\n";
        output << "                $isHeader = false;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $Value_Array, Round( $inFile_Line[1] * 1000 ));" << "\n";
        output << "            $Temp = Abs( Round( $inFile_Line[1] * 1000 ));" << "\n";
        output << "            if( $Temp > $Max ) $Max = $Temp;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Make Bin Array =====================-->" << "\n";
        output << "" << "\n";
        output << "        $Data_Array = Array();" << "\n";
        output << "        $Bin = $Bin_Size * 1000;" << "\n";
        output << "" << "\n";
        output << "        $Data_Array[ '<0' ] = Array();" << "\n";
        output << "        $Data_Array[ '>0' ] = Array();" << "\n";
        output << "" << "\n";
        output << "        $xMax = 0;" << "\n";
        output << "        $Max = (( $Max + $Bin ) - (( $Max + $Bin ) % $Bin ));" << "\n";
        output << "" << "\n";
        output << "        For( $i = $Max ; $i > 0; $i -= $Bin )" << "\n";
        output << "        {" << "\n";
        output << "            $Range = Number_Format((( $i - $Bin ) / 1000.0 ), 3, '.', '' ).'~'.Number_Format(( $i / 1000.0 ), 3, '.', '' );" << "\n";
        output << "            $Data_Array[ '<0' ][ $Range ] = 0;" << "\n";
        output << "            $Data_Array[ '>0' ][ $Range ] = 0;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Value_Array as $Value )" << "\n";
        output << "        {" << "\n";
        output << "            $Temp = Abs( $Value );" << "\n";
        output << "            $Temp = (( $Temp + $Bin ) - (( $Temp + $Bin ) % $Bin ));" << "\n";
        output << "            $Range = Number_Format((( $Temp - $Bin ) / 1000.0 ), 3, '.', '' ).'~'.Number_Format(( $Temp / 1000.0 ), 3, '.', '' );" << "\n";
        output << "            if( $Value < 0 ) $Data_Array[ '<0' ][ $Range ] -= 1;" << "\n";
        output << "            if( $Value > 0 ) $Data_Array[ '>0' ][ $Range ] += 1;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Data_Array as $Key => $Data )" << "\n";
        output << "        {" << "\n";
        output << "            Foreach( $Data as $Label => $Value )" << "\n";
        output << "            {" << "\n";
        output << "                if( $xMax < Abs( $Value )) $xMax = Abs( $Value );" << "\n";
        output << "            }" << "\n";
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
        output << "                    .margin({left: 80})" << "\n";
        output << "                    .groupSpacing(0)" << "\n";
        output << "                    .showControls(true)" << "\n";
        output << "                    .showValues(true)" << "\n";
        output << "                    .stacked(true)" << "\n";
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
