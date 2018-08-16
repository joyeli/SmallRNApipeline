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

    static void output_lendist_isomirs(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::ofstream output( output_name + sample_name + "-isomiRs.tsv" );
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

    static void output_lendist_trimmed(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::ofstream output( output_name + sample_name + "-trimmed.tsv" );
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

    static void output_lendist(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::vector< CountingTableType >& anno_table_trim, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name
            )
    {
        std::ofstream output( output_name + sample_name + ".tsv" );
        output << "Annotation\tA_Tail\tC_Tail\tG_Tail\tT_Tail\tOther_Tail\tGM";

        std::vector< std::string > split;
        std::set< std::string > anno_mark_set;
        std::set< std::string > anno_names;

        std::vector< CountingTableType > annos = std::vector< CountingTableType >( anno_table_tail.size() );
        std::vector< CountingTableType > whole = std::vector< CountingTableType >( anno_table_tail.size() );
        std::vector< CountingTableType > trims = std::vector< CountingTableType >( anno_table_tail.size() );

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            anno_names.emplace( split[0] );

            if( anno_mark.find( anno ) != anno_mark.end() )
                if( anno_mark[ anno ].at( anno_mark[ anno ].length() -1 ) == '!' )
                    anno_mark_set.emplace( split[0] );

            for( auto& len : ano_len_idx.second )
            {
                for( std::size_t i = 0; i < anno_table_tail.size(); ++i )
                {
                    if( annos[i][ split[0] ].find( len ) == annos[i][ split[0] ].end() ) annos[i][ split[0] ][ len ] = 0.0;
                    if( whole[i][ "allmir" ].find( len ) == whole[i][ "allmir" ].end() ) whole[i][ "allmir" ][ len ] = 0.0;
                    if( trims[i][ "allmir" ].find( len ) == trims[i][ "allmir" ].end() ) trims[i][ "allmir" ][ len ] = 0.0;

                    annos[i][ split[0] ][ len ] += anno_table_tail[i][ anno ][ len ];
                    whole[i][ "allmir" ][ len ] += anno_table_tail[i][ anno ][ len ];
                    trims[i][ "allmir" ][ len ] += anno_table_trim[i][ anno ][ len ];
                }
            }
        }

        for( auto& len : ano_len_idx.second )
        {
            output << "\n" << "All:" << len;

            for( auto& anno_table : whole )
            {
                if( anno_table.find( "allmir" ) != anno_table.end() )
                {
                    if( anno_table[ "allmir" ].find( len ) != anno_table[ "allmir" ].end() )
                        output << "\t" << std::setprecision( 0 ) << std::fixed << anno_table[ "allmir" ][ len ];
                    else output << "\t0";
                }
                else output << "\t0";
            }
        }

        for( auto& len : ano_len_idx.second )
        {
            output << "\n" << "Alltrim:" << len;

            for( auto& anno_table : trims )
            {
                if( anno_table.find( "allmir" ) != anno_table.end() )
                {
                    if( anno_table[ "allmir" ].find( len ) != anno_table[ "allmir" ].end() )
                        output << "\t" << std::setprecision( 0 ) << std::fixed << anno_table[ "allmir" ][ len ];
                    else output << "\t0";
                }
                else output << "\t0";
            }
        }

        for( auto& anno : anno_names )
        {
            for( auto& len : ano_len_idx.second )
            {
                output
                    << "\n" << anno << ( anno_mark_set.find( anno ) != anno_mark_set.end() ? "!" : "" )
                    << ":" << len;

                for( auto& anno_table : annos )
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

    static void output_lendist_visualization( const std::string& output_name, const bool& isSeed, const bool is_biotype = false )
    {
        std::ofstream output( output_name + "index.php" );
        bool trimming = false;

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "" << "\n";
        output << "        $ForceY = $_POST['ForceY'];" << "\n";
        output << "        $IsomiRs = $_POST['IsomiRs'];" << "\n";
        output << "        $Trimmed = $_POST['Trimmed'];" << "\n";
        output << "        $PPMFilter = $_POST['PPMFilter'];" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = $_GET['TSV_File'];" << "\n";
        output << "        $annotation_select = $_GET['annotation_select'];" << "\n";
        output << "        $isAbundant = $IsomiRs == 'No' ? 'AllAnnotations' : $_POST['isAbundant'];" << "\n";
        output << "" << "\n";
        output << "        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        if( $_GET['annotation_select'] != '' )" << "\n";
        output << "        {" << "\n";
        output << "            $IsomiRs = 'Yes';" << "\n";
        output << "            $isAbundant = 'AllAnnotations';" << "\n";
        output << "        }" << "\n";
        output << "        else $annotation_select = $_POST['annotation_select'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js ></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";

        if( is_biotype )
        {
            output << "#<--================== IsTrimmed ====================-->" << "\n";
            output << "" << "\n";
            output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
            output << "" << "\n";
            output << "        echo '<select name=Trimmed onchange=this.form.submit();>';" << "\n";
            output << "        echo '<option '; if($Trimmed=='') echo 'selected'; echo '>Is Trimming?</option>';" << "\n";
            output << "" << "\n";
            output << "        $miR_List = array('Yes', 'No');" << "\n";
            output << "" << "\n";
            output << "        For( $i = 0; $i < Count( $miR_List ); ++$i )" << "\n";
            output << "        {" << "\n";
            output << "            echo '<option value='.$miR_List[$i].' ';" << "\n";
            output << "" << "\n";
            output << "            if( $Trimmed == $miR_List[$i] )" << "\n";
            output << "                echo 'selected ';" << "\n";
            output << "" << "\n";
            output << "            echo '>' . $miR_List[$i] . '</option>';" << "\n";
            output << "        }" << "\n";
            output << "" << "\n";
            output << "        echo \"</select>" << "\n";
            output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
            output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
            output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
            output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
            output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
            output << "            </form>\";" << "\n";
            output << "" << "\n";
        }
        else output << "        $Trimmed = 'No';" << "\n";

        if( !isSeed && !is_biotype )
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
            output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
            output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
            output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
            output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
            output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
            output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
            output << "            </form>\";" << "\n";
        }
        else if( isSeed )
        {
            output << "        $IsomiRs = 'No';" << "\n";
        }
        else if( is_biotype )
        {
            output << "        $IsomiRs = 'Yes';" << "\n";
        }

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
        output << "        $TSV_List_Temp = array();" << "\n";
        output << "" << "\n";

        if( is_biotype )
        {
            output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
            output << "        {" << "\n";
            output << "            $TSV_Temp = Explode( $Trimmed == 'Yes' ? '-trimmed' : '-isomiRs', $TSV_List[$i] );" << "\n";
            output << "" << "\n";
            output << "            if( Count( $TSV_Temp ) > 1 )" << "\n";
            output << "                Array_Push( $TSV_List_Temp, $TSV_Temp[0] );" << "\n";
            output << "        }" << "\n";
        }
        else
        {
            output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
            output << "        {" << "\n";
            output << "            $TSV_Temp = Explode( '-isomiRs', $TSV_List[$i] );" << "\n";
            output << "" << "\n";
            output << "            if( $IsomiRs == 'Yes' )" << "\n";
            output << "            {" << "\n";
            output << "                if( Count( $TSV_Temp ) > 1 )" << "\n";
            output << "                    Array_Push( $TSV_List_Temp, $TSV_Temp[0] );" << "\n";
            output << "            }" << "\n";
            output << "            else" << "\n";
            output << "            {" << "\n";
            output << "                if( Count( $TSV_Temp ) == 1 )" << "\n";
            output << "                    Array_Push( $TSV_List_Temp, Substr( $TSV_Temp[0], 0, Strlen( $TSV_Temp[0] ) - 4 ));" << "\n";
            output << "            }" << "\n";
            output << "        }" << "\n";
        }

        output << "" << "\n";
        output << "        $TSV_List = $TSV_List_Temp;" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";

        if( is_biotype )
        {
            output << "        $TSV_File_Temp = $TSV_File.( $Trimmed == 'Yes' ? '-trimmed.tsv' : '-isomiRs.tsv' );" << "\n";
        }
        else
        {
            output << "        $TSV_File_Temp = $TSV_File.( $IsomiRs == 'Yes' ? '-isomiRs.tsv' : '.tsv' );" << "\n";
        }

        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";

        if( !is_biotype && !isSeed )
        {
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
            output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
            output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
            output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
            output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
            output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
            output << "            </form>\";" << "\n";
            output << "" << "\n";
        }
        else
        {
            output << "        $isAbundant = 'AllAnnotations';" << "\n";
            output << "" << "\n";
        }

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
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( $TSV_File_Temp );" << "\n";
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
        output << "            if( $isAbundant == 'MostAbundant' && Count( $Annotation ) != 2 ) continue;" << "\n";
        output << "            if( !Array_Key_Exists( $anno[0], $Anno_Array )) $Anno_Array[ $anno[0] ] = 0;" << "\n";
        output << "" << "\n";
        output << "            for( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n";
        output << "                $Anno_Array[ $anno[0] ] += $inFile_Line[$k];" << "\n";
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
        output << "" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $is_any_selected = false;" << "\n";
        output << "        $Total_PPM = 0.0;" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "        {" << "\n";
        output << "            if( $PPMFilter != 'PPM Filter' && $PPMFilter != '' && $PPMFilter > $ppm )" << "\n";
        output << "                Unset( $Anno_Array[ $anno ] );" << "\n";
        output << "            else if( $Total_PPM < $ppm ) $Total_PPM = $ppm;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=annotation_select onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($annotation_select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$anno.' ';" << "\n";
        output << "" << "\n";
        output << "            if( $annotation_select == $anno ) " << "\n";
        output << "            {" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "                $is_any_selected = true;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$anno.' ('.$ppm.'ppm)</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        if( !$is_any_selected ) $Annotation_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='ForceY' value='$ForceY' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
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
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='Trimmed' value='$Trimmed' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='PPMFilter' value='$PPMFilter' />" << "\n";
        output << "            <input type='hidden' name='annotation_select' value='$annotation_select' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        
        if( !is_biotype && !isSeed )
        {
            output << "#<!--================== GoBack ====================-->" << "\n";
            output << "" << "\n";
            output << "        $annos = Explode( '-', $annotation_select );" << "\n";
            output << "        $anno  = $annos[0];" << "\n";
            output << "" << "\n";
            output << "        For( $i = 1; $i < Count( $annos ) -1; $i++ ) $anno = $anno.'-'.$annos[ $i ];" << "\n";
            output << "" << "\n";
            output << "        $annos = Explode( '_', $anno );" << "\n";
            output << "        $anno  = $annos[0];" << "\n";
            output << "" << "\n";
            output << "        For( $i = 1; $i < Count( $annos ) -1; $i++ ) $anno = $anno.'_'.$annos[ $i ];" << "\n";
            output << "" << "\n";
            output << "        if( $TSV_File != '' && $annotation_select != '' )" << "\n";
            output << "            echo \"<a target='_blank' href='../SqAlign/index.php?TSV_File=$TSV_File.tsv&Annotation_Select=$anno' >" << "\n";
            output << "                <input type='submit' value='Goto $anno' />" << "\n";
            output << "                </a>\";" << "\n";
            output << "" << "\n";
        }

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
