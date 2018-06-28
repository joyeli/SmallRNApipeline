#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBarplot
{

  public:

    GeneTypeAnalyzerBarplot()
    {}

    static void output_barplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& biotype,
            const std::string& token
            )
    {
        std::ofstream output;
        std::vector< std::pair< std::string, double >> total_vec;
        std::map< std::string, std::map< std::string, std::vector< double >>> mirlen_map;

        make_table( bed_samples, ano_len_idx, anno_table_tail, token, total_vec, mirlen_map );

        boost::filesystem::create_directory( boost::filesystem::path( output_name + "Loading/" ));
        sort_total_vec( total_vec );

        for( auto& len : ano_len_idx.second )
        {
            output.open( output_name + "Loading/" + token + "_" + std::to_string( len ) + ".tsv" );
            output_table( output, std::to_string( len ), bed_samples, total_vec, mirlen_map );
            output.close();
        }

        boost::filesystem::create_directory( boost::filesystem::path( output_name + "Length/" ));
        merge_seed( bed_samples, total_vec, mirlen_map );
        sort_total_vec( total_vec );

        for( auto& len : ano_len_idx.second )
        {
            output.open( output_name + "Length/" + token + "_" + std::to_string( len ) + ".tsv" );
            output_table( output, std::to_string( len ), bed_samples, total_vec, mirlen_map );
            output.close();
        }

        if( biotype.substr( 0, 5 ) == "miRNA" || biotype == "mirtron" )
        {
            boost::filesystem::create_directory( boost::filesystem::path( output_name + "Arm/" ));
            merge_length( bed_samples, total_vec, mirlen_map );
            sort_total_vec( total_vec );

            output.open( output_name + "Arm/" + token + "_5p" + ".tsv" );
            output_table( output, "5p", bed_samples, total_vec, mirlen_map );
            output.close();

            output.open( output_name + "Arm/" + token + "_3p" + ".tsv" );
            output_table( output, "3p", bed_samples, total_vec, mirlen_map );
            output.close();
        }
    }

    static void make_table(
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const std::string& token,
            std::vector< std::pair< std::string, double >>& total_vec,
            std::map< std::string, std::map< std::string, std::vector< double >>>& mirlen_map
            )
    {
        double gm = 0.0;
        double pm = 0.0;
        std::map< std::size_t, std::vector< std::pair< double, double >>> len_map;

        for( auto& anno : ano_len_idx.first )
        {
            gm = 0.0;
            pm = 0.0;
            len_map.clear();

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                    for( auto& len : ano_len_idx.second )
                    {
                        if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                            gm += anno_table_tail[ smp ][5][ anno ][ len ];

                        if( len_map.find( len ) == len_map.end() )
                            len_map[ len ] = std::vector< std::pair< double, double >>( bed_samples.size(), std::make_pair( 0.0, 0.0 ));

                        len_map[ len ][ smp ].first += anno_table_tail[ smp ][5][ anno ][ len ];
                    }

                for( std::size_t i = 0; i < 5; i++ )
                {
                    if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                pm += anno_table_tail[ smp ][i][ anno ][ len ];

                            if( len_map.find( len ) == len_map.end() )
                                len_map[ len ] = std::vector< std::pair< double, double >>( bed_samples.size(), std::make_pair( 0.0, 0.0 ));

                            len_map[ len ][ smp ].second += anno_table_tail[ smp ][i][ anno ][ len ];
                        }
                }

                for( auto& len : len_map )
                {
                    if( mirlen_map[ anno ].find( std::to_string( len.first )) == mirlen_map[ anno ].end() )
                        mirlen_map[ anno ][ std::to_string( len.first )] = std::vector< double >( bed_samples.size(), 0.0 );
                    mirlen_map[ anno ][ std::to_string( len.first )][ smp ] = 
                        ( token == "GMPM"
                        ? len.second[ smp ].first + len.second[ smp ].second
                        : ( token == "GM"
                            ? len.second[ smp ].first
                            : ( token == "PM"
                                ? len.second[ smp ].second
                                : (( len.second[ smp ].first + len.second[ smp ].second ) < 1
                                    ? 0
                                    : ( len.second[ smp ].second * 100 / ( len.second[ smp ].first + len.second[ smp ].second ))))));
                }
            }

            total_vec.emplace_back( anno, gm + pm );
        }
    }

    static void merge_seed(
            const std::vector< BedSampleType >& bed_samples,
            std::vector< std::pair< std::string, double >>& total_vec,
            std::map< std::string, std::map< std::string, std::vector< double >>>& mirlen_map
            )
    {
        std::string name;
        std::map< std::string, double > total_tmp;
        std::map< std::string, std::map< std::string, std::vector< double >>> mirlen_tmp;

        for( auto& anno : total_vec )
        {
            name = anno.first.substr( 0, anno.first.length() -8 );

            if( total_tmp.find( name ) == total_tmp.end() )
                total_tmp[ name ] = 0.0;

            total_tmp[ name ] += anno.second;

            for( auto& len : mirlen_map[ anno.first ] )
            {
                if( mirlen_tmp[ name ].find( len.first ) == mirlen_tmp[ name ].end() )
                    mirlen_tmp[ name ][ len.first ] = std::vector< double >( bed_samples.size(), 0.0 );
                
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                    mirlen_tmp[ name ][ len.first ][ smp ] += mirlen_map[ anno.first ][ len.first ][ smp ];
            }
        }

        total_vec.clear();
        mirlen_map = mirlen_tmp;

        for( auto& anno : total_tmp )
            total_vec.emplace_back( anno );
    }

    static void merge_length(
            const std::vector< BedSampleType >& bed_samples,
            std::vector< std::pair< std::string, double >>& total_vec,
            std::map< std::string, std::map< std::string, std::vector< double >>>& mirlen_map
            )
    {
        std::string arm, name;
        std::map< std::string, double > total_tmp;
        std::map< std::string, std::map< std::string, std::vector< double >>> mirlen_tmp;

        for( auto& anno : total_vec )
        {
            arm  = anno.first.substr( anno.first.length() -2, 2 );
            name = anno.first.substr( 0, anno.first.length() -3 );

            if( total_tmp.find( name ) == total_tmp.end() )
                total_tmp[ name ] = 0.0;

            total_tmp[ name ] += anno.second;

            if( mirlen_tmp[ name ].find( arm ) == mirlen_tmp[ name ].end() )
                mirlen_tmp[ name ][ arm ] = std::vector< double >( bed_samples.size(), 0.0 );
                
            for( auto& len : mirlen_map[ anno.first ] )
            {
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                    mirlen_tmp[ name ][ arm ][ smp ] += mirlen_map[ anno.first ][ len.first ][ smp ];
            }
        }

        total_vec.clear();
        mirlen_map = mirlen_tmp;

        for( auto& anno : total_tmp )
            total_vec.emplace_back( anno );
    }

    static void sort_total_vec( std::vector< std::pair< std::string, double >>& total_vec )
    {
        std::sort( total_vec.begin(), total_vec.end(),
        []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
        {
            if( a.second == b.second )
                return a.first > b.first;
            else
                return a.second > b.second;
        });
    }

    static void output_table(
            std::ofstream& output,
            const std::string tag,
            const std::vector< BedSampleType >& bed_samples,
            const std::vector< std::pair< std::string, double >>& total_vec,
            std::map< std::string, std::map< std::string, std::vector< double >>>& mirlen_map
            )
    {
        output << "Annotation";
        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : total_vec )
        {
            output << "\n" << anno.first << std::setprecision( 0 ) << std::fixed;
            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( mirlen_map[ anno.first ].find( tag ) != mirlen_map[ anno.first ].end() )
                    output << "\t" << mirlen_map[ anno.first ][ tag ][ smp ];
                else
                    output << "\t" << 0;
        }
        output << "\n";
    }

    static void output_barplot_visualization(
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
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $MaxHight = $_POST['MaxHight'];" << "\n";
        output << "        $Top_miRNA = $_POST['Top_miRNA'];" << "\n";
        output << "        $Min_Length = $_POST['Min_Length'];" << "\n";
        output << "        $Max_Length = $_POST['Max_Length'];" << "\n";
        output << "        $barPlotType = $_POST['barPlotType'];" << "\n";
        output << "        $Select_Type = $_POST['Select_Type'];" << "\n";
        output << "        $Single_Anno = $_POST['Single_Anno'];" << "\n";
        output << "        $Filter_Option = $_POST['Filter_Option'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/d3.min.js></script>';" << "\n";
        output << "        echo '<link href=http://192.168.1.11:6680/ness/www/ForAgoSorting/lib/svg0331.css rel=stylesheet type=text/css>';" << "\n";
        output << "        echo '<style> .x.axis path { display: none; }</style>';" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "#<!--================= Filter ===================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Filter onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($Filter=='') echo 'selected'; echo '>noFilter</option>';" << "\n";
        output << "" << "\n";
        output << "        $Filters = Shell_Exec( 'ls ./ | grep -v php' );" << "\n";
        output << "        $Filter_List = Explode( \"\\n\", $Filters );" << "\n";
        output << "        $Filter_Size  = Count( $Filter_List )-1;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Filter_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Filter_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Filter == $Filter_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $Filter_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--============== Filter Option =================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == 'Length' || $Filter == 'Loading' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=Filter_Option size=5 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $Filter_Option == '' || $Filter_Option == 'LenDiff' || $Filter_Option == '%Change' )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Filter == 'Length'  ) echo 'LenDiff';" << "\n";
        output << "                if( $Filter == 'Loading' ) echo '%Change';" << "\n";
        output << "            }" << "\n";
        output << "            else echo $Filter_Option;" << "\n";
        output << "" << "\n";
        output << "            echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=GMPM onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($GMPM=='') echo 'selected'; echo '>GMPM</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPM_List = array('GM', 'PM', 'Tailing');" << "\n";
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
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== barPlotType ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=barPlotType onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $barPlotType_List = array('ppm', '100%');" << "\n";
        output << "        $barPlotType_Size = Count( $barPlotType_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $barPlotType_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$barPlotType_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $barPlotType == $barPlotType_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $barPlotType_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Max NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $barPlotType != '100%' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=MaxHight size=8 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $MaxHight=='' )" << "\n";
        output << "                echo 'MaxHight';" << "\n";
        output << "            else" << "\n";
        output << "                echo $MaxHight;" << "\n";
        output << "" << "\n";
        output << "            echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Select Type ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Select_Type onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $Select_Type_List = array('SelectType', 'SingleAnno', 'TopRanking');" << "\n";
        output << "        $Select_Type_Size = Count( $Select_Type_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $Select_Type_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Select_Type_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Select_Type == $Select_Type_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $Select_Type_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Top miRNA Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Select_Type == 'TopRanking' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=Top_miRNA size=4 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $Top_miRNA == '' )" << "\n";
        output << "            {" << "\n";
        output << "                echo 'miRNA#';" << "\n";
        output << "            }" << "\n";
        output << "            else echo $Top_miRNA;" << "\n";
        output << "" << "\n";
        output << "            echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Getting Filter ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Filtered_miRNAs = Array();" << "\n";
        output << "        $Filter_File = File_get_contents( '../Difference/'.$Filter.'Difference_'.$GMPM.'.tsv' );" << "\n";
        output << "        $Filter_Line = Explode( \"\\n\", $Filter_File );" << "\n";
        output << "        $count = -1;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 1; $i < Count( $Filter_Line )-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $Line_Splits = Explode( \"\\t\", $Filter_Line[$i] );" << "\n";
        output << "" << "\n";
        output << "            if( $Filter == 'Arm' )" << "\n";
        output << "                if( $Line_Splits[2] == 'Y' )" << "\n";
        output << "                {" << "\n";
        output << "                    $count++;" << "\n";
        output << "                    if( $Top_miRNA != '' && $count >= $Top_miRNA+1 ) break;" << "\n";
        output << "                    Array_Push( $Filtered_miRNAs, $Line_Splits[0] );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "            if( $Filter == 'Length' && $Filter_Option != 'LenDiff' )" << "\n";
        output << "                if( $Line_Splits[2] >= $Filter_Option )" << "\n";
        output << "                {" << "\n";
        output << "                    $count++;" << "\n";
        output << "                    if( $Top_miRNA != '' && $count >= $Top_miRNA+1 ) break;" << "\n";
        output << "                    Array_Push( $Filtered_miRNAs, $Line_Splits[0] );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "            if( $Filter == 'Loading' && $Filter_Option != '%Change' )" << "\n";
        output << "                if(( $Line_Splits[2] > 0 && $Line_Splits[2] >= $Filter_Option ) || ( $Line_Splits[2] < 0 && $Line_Splits[2] <= $Filter_Option *-1 ))" << "\n";
        output << "                {" << "\n";
        output << "                    $count++;" << "\n";
        output << "                    if( $Top_miRNA != '' && $count >= $Top_miRNA+1 ) break;" << "\n";
        output << "                    Array_Push( $Filtered_miRNAs, $Line_Splits[0] );" << "\n";
        output << "                }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== Single Annotation ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Select_Type == 'SingleAnno' )" << "\n";
        output << "        {" << "\n";
        output << "            $Single_Anno_List = array();" << "\n";
        output << "" << "\n";
        output << "            $inFile = File_get_contents(( $Filter == 'noFilter' ? 'Loading' : $Filter ).'/'.$GMPM.'_'.( $Filter == 'Arm' ? '5p' : '20' ).'.tsv' );" << "\n";
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "            For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "" << "\n";
        output << "                if( Count( $Filtered_miRNAs ) != 0 && !in_array( $inFile_Line[0], $Filtered_miRNAs ))" << "\n";
        output << "                    continue;" << "\n";
        output << "" << "\n";
        output << "                Array_Push( $Single_Anno_List, $inFile_Line[0] );" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '<select name=Single_Anno onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "            Rsort( $Single_Anno_List );" << "\n";
        output << "            $Single_Anno_Size = Count( $Single_Anno_List );" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < $Single_Anno_Size; ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$Single_Anno_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $Single_Anno == $Single_Anno_List[$i] )" << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "                echo '>' . $Single_Anno_List[$i] . '</option>';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>\";" << "\n";
        output << "        }" << "\n";
        output << "        else $Single_Anno = '';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MinLength input =====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Filter != 'Arm' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=Min_Length size=3 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $Min_Length == '' )" << "\n";
        output << "            {" << "\n";
        output << "                echo 'minLen';" << "\n";
        output << "            }" << "\n";
        output << "            else echo $Min_Length;" << "\n";
        output << "" << "\n";
        output << "            echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MaxLength input =====================-->" << "\n";
        output << "" << "\n";
        output << "            echo '<input type=text onchange=this.form.submit(); name=Max_Length size=3 value=';" << "\n";
        output << "" << "\n";
        output << "            if( $Max_Length == '' )" << "\n";
        output << "            {" << "\n";
        output << "                echo 'maxLen';" << "\n";
        output << "            }" << "\n";
        output << "            else echo $Max_Length;" << "\n";
        output << "" << "\n";
        output << "            echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"<input type='submit' value='Submit' />\";" << "\n";
        output << "        echo '</form><br/>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Multi TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = Array();" << "\n";
        output << "" << "\n";
        output << "        if( $Filter != 'Arm' )" << "\n";
        output << "        {" << "\n";
        output << "            For( $i = $Min_Length; $i <= $Max_Length; ++$i )" << "\n";
        output << "                Array_Push( $TSV_File, ( $Filter == 'noFilter' ? 'Loading' : $Filter ).'/'.$GMPM.'_'.$i.'.tsv' );" << "\n";
        output << "        }" << "\n";
        output << "        else" << "\n";
        output << "        {" << "\n";
        output << "            Array_Push( $TSV_File, $Filter.'/'.$GMPM.'_5p.tsv' );" << "\n";
        output << "            Array_Push( $TSV_File, $Filter.'/'.$GMPM.'_3p.tsv' );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--================== barPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        //========== svg view var set ==========" << "\n";
        output << "" << "\n";
        output << "        $width = ( $Top_miRNA == '' ? 1 : $Top_miRNA ) * 180;" << "\n";
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
        output << "        if( $barPlotType == 'ppm' )" << "\n";
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
        output << "                        if( $Single_Anno != '' && $inFile_Line[0] != $Single_Anno ) continue;" << "\n";
        output << "                        if( Count( $Filtered_miRNAs ) != 0 && !in_array( $inFile_Line[0], $Filtered_miRNAs ))" << "\n";
        output << "                            continue;" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( $inFile_Line[$k] >= $MaxValue )" << "\n";
        output << "                                $MaxValue = $inFile_Line[$k];" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "            else $MaxValue = $MaxHight;" << "\n";
        output << "        }" << "\n";
        output << "        else $MaxValue = 100;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_File ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = File_get_contents( $TSV_File[$i] );" << "\n";
        output << "            $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "            $Temp = Tempnam( '/tmp', $TSV_File[$i] );" << "\n";
        output << "            $Ftemp = fopen( $Temp, 'w' );" << "\n";
        output << "            $count = -1;" << "\n";
        output << "" << "\n";
        output << "            if( $barPlotType == 'ppm' )" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $count++;" << "\n";
        output << "                    if( $Single_Anno == '' && $count >= $Top_miRNA+1 )" << "\n";
        output << "                        break;" << "\n";
        output << "" << "\n";
        output << "                    if( $j == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n";
        output << "                        continue;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "" << "\n";
        output << "                    if( $Single_Anno != '' && $inFile_Line[0] != $Single_Anno ) continue;" << "\n";
        output << "                    if( Count( $Filtered_miRNAs ) != 0 && !in_array( $inFile_Line[0], $Filtered_miRNAs ))" << "\n";
        output << "                    {" << "\n";
        output << "                        $count--;" << "\n";
        output << "                        continue;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    $count++;" << "\n";
        output << "                    if( $Single_Anno == '' && $count >= $Top_miRNA+1 )" << "\n";
        output << "                        break;" << "\n";
        output << "" << "\n";
        output << "                    if( $j == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n";
        output << "                        continue;" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "" << "\n";
        output << "                        if( $Single_Anno != '' && $inFile_Line[0] != $Single_Anno ) continue;" << "\n";
        output << "                        if( Count( $Filtered_miRNAs ) != 0 && !in_array( $inFile_Line[0], $Filtered_miRNAs ))" << "\n";
        output << "                        {" << "\n";
        output << "                            $count--;" << "\n";
        output << "                            continue;" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        $Total_Value = 0;" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                            $Total_Value += $inFile_Line[ $k ];" << "\n";
        output << "" << "\n";
        output << "                        fwrite( $Ftemp, $inFile_Line[0] );" << "\n";
        output << "" << "\n";
        output << "                        For( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                            fwrite( $Ftemp, \"\\t\".Number_format( $inFile_Line[ $k ]*100/$Total_Value, 0 ));" << "\n";
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
        output << "                    .call(xAxis)" << "\n";
        output << "                    .selectAll('text')" << "\n";
        output << "                    .style('font-size', function(d){" << "\n";
        output << "                        return 75 - ( d.length > 15 ? (( d.length - 15 ) * 1.25 ) : 0 ) + '%';" << "\n";
        output << "                    });" << "\n";
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
        output << "                    .data(data).enter()" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', 18)" << "\n";
        output << "                    .attr('y1', 0)" << "\n";
        output << "                    .attr('x2', 18)" << "\n";
        output << "                    .attr('y2', height)" << "\n";
        output << "                    .attr('style', 'stroke:rgb(192,192,192);stroke-width:1')" << "\n";
        output << "                    .attr('transform', function(d) { return 'translate(' + (x0(d.Annotation) + d.values.length * x1.rangeBand()) + ',0)'; });" << "\n";
        output << "            " << "\n";
        output << "                    barchat = svg$i.selectAll('.barchat')" << "\n";
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
