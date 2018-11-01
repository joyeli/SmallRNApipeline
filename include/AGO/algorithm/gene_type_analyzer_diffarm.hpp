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
                diffarm[ anno_arm.first ] = std::vector< std::vector< double >>( 3, std::vector< double >( 18, 0.0 ));

            diffarm[ anno_arm.first ][ anno_arm.second == "5p" ? 0 : 1 ] = 
            {
                std::stod( split[1] ),
                0.0, // GM
                std::stod( split[1] ) * std::stod( split[3] ),
                std::stod( split[3] ),
                std::stod( split[4] ),
                std::stod( split[5] ),
                std::stod( split[6] ),
                std::stod( split[7] ),
                std::stod( split[8] ),
                0.0, // RatioGMPM
                0.0, // RatioGM
                0.0, // RatioPM
                0.0, // RatioTailing
                0.0, // RatioA
                0.0, // RatioC
                0.0, // RatioG
                0.0, // RatioU
                0.0  // RatioO
            };

            diffarm[ anno_arm.first ][ anno_arm.second == "5p" ? 0 : 1 ][1]
                = diffarm[ anno_arm.first ][ anno_arm.second == "5p" ? 0 : 1 ][0]
                - diffarm[ anno_arm.first ][ anno_arm.second == "5p" ? 0 : 1 ][2]
                ;
        }

        infile.close();
        return diffarm;
    }

    static void get_diff( DiffArmType& diffarm )
    {
        double tmp = 0.0;

        for( auto& anno : diffarm )
        {
            auto& p5 = anno.second[0];
            auto& p3 = anno.second[1];
            auto& p2 = anno.second[2];

            p2[0] = p5[0] + p3[0];
            p2[1] = p5[1] + p3[1];
            p2[2] = p5[2] + p3[2];
            p2[3] = p2[0] == 0 ? 0 : ( p2[2] / p2[0] );

            for( std::size_t i = 9; i < anno.second[2].size(); ++i )
            {
                tmp = p5[ i-9 ] + p3[ i-9 ];

                p5[i] = tmp == 0 ? 0 : ( p5[ i-9 ] / tmp );
                p3[i] = tmp == 0 ? 0 : ( p3[ i-9 ] / tmp );
            }

            for( std::size_t i = 4; i < 9; ++i )
                p2[i+1] = p5[i] - p3[i];
        }
    }

    static void out_diff( DiffArmType& diffarm, const std::string& out_file )
    {
        std::ofstream output( out_file );
        output << "Anno\t"
            << "GMPM\tGM\tPM\tTailing\tDiff_Tailing\tDiff_A\tDiff_C\tDiff_G\tDiff_U\tDiff_Other\t"
            << "GMPM\tGM\tPM\tTailing\tAtail\tCtail\tGtail\tUtail\tOther\t"
            << "Ratio_GMPM\tRatio_GM\tRatio_PM\tRatio_Tailing\tRatio_Atail\tRatio_Ctail\tRatio_Gtail\tRatio_Utail\tRatio_Other\t"
            << "GMPM\tGM\tPM\tTailing\tAtail\tCtail\tGtail\tUtail\tOther\t"
            << "Ratio_GMPM\tRatio_GM\tRatio_PM\tRatio_Tailing\tRatio_Atail\tRatio_Ctail\tRatio_Gtail\tRatio_Utail\tRatio_Other\n"
            ;

        for( auto& anno : diffarm )
        {
            output << anno.first;

            for( std::size_t i = 0; i < 10                   ; ++i ) output << "\t" << anno.second[2][i];
            for( std::size_t i = 0; i < anno.second[0].size(); ++i ) output << "\t" << anno.second[0][i];
            for( std::size_t i = 0; i < anno.second[1].size(); ++i ) output << "\t" << anno.second[1][i];

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
        output << "        //  0 Anno" << "\n";
        output << "        //  1 GMPM" << "\n";
        output << "        //  2 Tailing" << "\n";
        output << "        //  3 Diff_Tailing" << "\n";
        output << "        //  4 Diff_A" << "\n";
        output << "        //  5 Diff_C" << "\n";
        output << "        //  6 Diff_G" << "\n";
        output << "        //  7 Diff_U" << "\n";
        output << "        //  8 Diff_Other" << "\n";
        output << "        //  9 GMPM" << "\n";
        output << "        // 10 Tailing" << "\n";
        output << "        // 11 Atail" << "\n";
        output << "        // 12 Ctail" << "\n";
        output << "        // 13 Gtail" << "\n";
        output << "        // 14 Utail" << "\n";
        output << "        // 15 Other" << "\n";
        output << "        // 16 Ratio5p" << "\n";
        output << "        // 17 GMPM" << "\n";
        output << "        // 18 Tailing" << "\n";
        output << "        // 19 Atail" << "\n";
        output << "        // 20 Ctail" << "\n";
        output << "        // 21 Gtail" << "\n";
        output << "        // 22 Utail" << "\n";
        output << "        // 23 Other" << "\n";
        output << "        // 24 Ratio3p" << "\n";
        output << "" << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "" << "\n";
        output << "        $isPair = $_POST['isPair'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $Filter5p = $_POST['Filter5p'];" << "\n";
        output << "        $Filter3p = $_POST['Filter3p'];" << "\n";
        output << "        $DiffType = $_POST['DiffType'];" << "\n";
        output << "        $isTrimmed = $_POST['isTrimmed'];" << "\n";
        output << "        $FilterGMPM = $_POST['FilterGMPM'];" << "\n";
        output << "        $FilteRatio = $_POST['FilteRatio'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.6/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
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
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Filter5p' value='$Filter5p' />" << "\n";
        output << "            <input type='hidden' name='Filter3p' value='$Filter3p' />" << "\n";
        output << "            <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='hidden' name='FilterGMPM' value='$FilterGMPM' />" << "\n";
        output << "            <input type='hidden' name='FilteRatio' value='$FilteRatio' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== TSV Files ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $TSVs = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $TSV_List ); ++$i )" << "\n";
        output << "            if( $TSV_List[$i] != '' ) Array_Push( $TSVs, $TSV_List[$i] );" << "\n";
        output << "" << "\n";
        output << "        if( $isPair == 'Yes' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "            echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "            echo '<option '; if($TSV_File=='') echo 'selected'; echo '>TSV_File</option>';" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $TSVs ); ++$i )" << "\n";
        output << "                echo '<option value='.$TSVs[$i].( $TSV_File == $TSVs[$i] ? ' selected >' : ' >' ).$TSVs[$i].'</option>';" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "                <input type='hidden' name='Filter5p' value='$Filter5p' />" << "\n";
        output << "                <input type='hidden' name='Filter3p' value='$Filter3p' />" << "\n";
        output << "                <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "                <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "                <input type='hidden' name='FilterGMPM' value='$FilterGMPM' />" << "\n";
        output << "                <input type='hidden' name='FilteRatio' value='$FilteRatio' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "" << "\n";
        output << "            $TSVs = Array();" << "\n";
        output << "            if( $TSV_File != '' && $TSV_File != 'TSV_File' )" << "\n";
        output << "                Array_Push( $TSVs, $TSV_File );" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Get Header ====================-->" << "\n";
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
        output << "                For( $i = 3; $i <= 8; ++$i )" << "\n";
        output << "                    $Header[ $inFile_Line[$i] ] = $i;" << "\n";
        output << "" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "            else break;" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Header[ 'Ratios_5p3p' ] = -1;" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=DiffType onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($DiffType=='') echo 'selected'; echo 'value= >DiffType</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Header as $Type => $Index )" << "\n";
        output << "            echo '<option value='.$Type.' '.( $DiffType == $Type ? ' selected >' : ' >' ).$Type.'</option>';" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Filter5p' value='$Filter5p' />" << "\n";
        output << "            <input type='hidden' name='Filter3p' value='$Filter3p' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='hidden' name='FilterGMPM' value='$FilterGMPM' />" << "\n";
        output << "            <input type='hidden' name='FilteRatio' value='$FilteRatio' />" << "\n";
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
        output << "            echo '<option value='.$Trim_List[$i].( $isTrimmed == $Trim_List[$i] ? ' selected >±' : ' >±' ).$Trim_List[$i].'％</option>';" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Filter5p' value='$Filter5p' />" << "\n";
        output << "            <input type='hidden' name='Filter3p' value='$Filter3p' />" << "\n";
        output << "            <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "            <input type='hidden' name='FilterGMPM' value='$FilterGMPM' />" << "\n";
        output << "            <input type='hidden' name='FilteRatio' value='$FilteRatio' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Filtering =====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $FilterGMPM == '' ) $FilterGMPM = 'FilterGMPM';" << "\n";
        output << "        if( $FilteRatio == '' ) $FilteRatio = 'Ratio%';" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=FilterGMPM size=7 value='.$FilterGMPM.' onfocus=\"{this.value=\\'\\';}\">';" << "\n";
        output << "        echo '<input type=text name=FilteRatio size=3 value='.$FilteRatio.' onfocus=\"{this.value=\\'\\';}\">';" << "\n";
        output << "" << "\n";
        output << "        $List5p = Array();" << "\n";
        output << "        $List3p = Array();" << "\n";
        output << "" << "\n";
        output << "        $List5p[ 'Atail' ] = 11;" << "\n";
        output << "        $List5p[ 'Ctail' ] = 12;" << "\n";
        output << "        $List5p[ 'Gtail' ] = 13;" << "\n";
        output << "        $List5p[ 'Utail' ] = 14;" << "\n";
        output << "        $List5p[ 'Other' ] = 15;" << "\n";
        output << "" << "\n";
        output << "        $List3p[ 'Atail' ] = 19;" << "\n";
        output << "        $List3p[ 'Ctail' ] = 20;" << "\n";
        output << "        $List3p[ 'Gtail' ] = 21;" << "\n";
        output << "        $List3p[ 'Utail' ] = 22;" << "\n";
        output << "        $List3p[ 'Other' ] = 23;" << "\n";
        output << "" << "\n";
        output << "        if( $DiffType == 'Ratios_5p3p' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<select name=Filter5p>';" << "\n";
        output << "            echo '<option '; if($Filter5p=='') echo 'selected'; echo '>Filter5p</option>';" << "\n";
        output << "" << "\n";
        output << "            Foreach( $List5p as $Tail => $Value )" << "\n";
        output << "                echo '<option value='.$Tail.( $Filter5p == $Tail ? ' selected >' : ' >' ).$Tail.'</option>';" << "\n";
        output << "" << "\n";
        output << "            echo '</select> > ';" << "\n";
        output << "            echo '<select name=Filter3p>';" << "\n";
        output << "            echo '<option '; if($Filter3p=='') echo 'selected'; echo '>Filter3p</option>';" << "\n";
        output << "" << "\n";
        output << "            Foreach( $List3p as $Tail => $Value )" << "\n";
        output << "                echo '<option value='.$Tail.( $Filter3p == $Tail ? ' selected >' : ' >' ).$Tail.'</option>';" << "\n";
        output << "" << "\n";
        output << "            echo '</select>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='isPair' value='$isPair' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='DiffType' value='$DiffType' />" << "\n";
        output << "            <input type='hidden' name='isTrimmed' value='$isTrimmed' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "        $FilteRatio = $FilteRatio / 100;" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Read File =====================-->" << "\n";
        output << "       " << "\n";
        output << "        $yMin = 0;" << "\n";
        output << "        $yMax = 0;" << "\n";
        output << "        $Tail_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        if( $DiffType != '' && $DiffType != 'DiffType' )" << "\n";
        output << "        {" << "\n";
        output << "            For( $i = 0; $i < Count( $TSVs ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $Tail_Array [$i] = Array();" << "\n";
        output << "" << "\n";
        output << "                $isHeader = true;" << "\n";
        output << "                $inFile = new SplFileObject( $TSVs[$i] );" << "\n";
        output << "" << "\n";
        output << "                while( !$inFile->eof() )" << "\n";
        output << "                {" << "\n";
        output << "                    $inFile_Lines = $inFile->fgets();" << "\n";
        output << "" << "\n";
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
        output << "" << "\n";
        output << "                    if( $FilterGMPM != '' && $FilterGMPM != 'FilterGMPM' && (Float) $inFile_Line[1] <= $FilterGMPM ) continue;" << "\n";
        output << "                    if( $FilteRatio != '' && $FilteRatio != 'Ratio%'     && (Float) $inFile_Line[2] <= $FilteRatio ) continue;" << "\n";
        output << "" << "\n";
        output << "                    if( $isPair == 'Yes' )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Temp = Array();" << "\n";
        output << "                        $Name = $inFile_Line[0].'\\t'.Number_Format( $inFile_Line[1] ).'ppm\\t'.Number_Format( $inFile_Line[2] * 100 ).'%';" << "\n";
        output << "" << "\n";
        output << "                        if( $DiffType == 'Ratios_5p3p' )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( $Filter5p != $Filter3p && $Filter5p != 'Filter5p' && $Filter3p != 'Filter3p' )" << "\n";
        output << "                                if( $inFile_Line[ $List5p[ $Filter5p ]] <= $inFile_Line[ $List3p[ $Filter3p ]] ) continue;" << "\n";
        output << "" << "\n";
        output << "                            $Temp = Array(" << "\n";
        output << "                                $Name," << "\n";
        output << "                                $inFile_Line[ 16 ]," << "\n";
        output << "                                (-$inFile_Line[ 24 ])," << "\n";
        output << "                            );" << "\n";
        output << "" << "\n";
        output << "                            Array_Push( $Tail_Array[$i], $Temp );" << "\n";
        output << "                        }" << "\n";
        output << "                        else" << "\n";
        output << "                        {" << "\n";
        output << "                            $Temp = Array(" << "\n";
        output << "                                $Name," << "\n";
        output << "                                $inFile_Line[ $Header[ $DiffType ] +  7 ]," << "\n";
        output << "                                (-$inFile_Line[ $Header[ $DiffType ] + 15 ])," << "\n";
        output << "                            );" << "\n";
        output << "" << "\n";
        output << "                            Array_Push( $Tail_Array[$i], $Temp );" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "                    else if( $DiffType == 'Ratios_5p3p' )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( $Filter5p != $Filter3p && $Filter5p != 'Filter5p' && $Filter3p != 'Filter3p' )" << "\n";
        output << "                           if( $inFile_Line[ $List5p[ $Filter5p ]] <= $inFile_Line[ $List3p[ $Filter3p ]] ) continue;" << "\n";
        output << "" << "\n";
        output << "                        Array_Push( $Tail_Array[$i], $inFile_Line[ 16 ]);" << "\n";
        output << "                        Array_Push( $Tail_Array[$i], $inFile_Line[ 24 ]);" << "\n";
        output << "                    }" << "\n";
        output << "                    else Array_Push( $Tail_Array[$i], $inFile_Line[ $Header[ $DiffType ]] );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        if( $isPair == 'Yes' )" << "\n";
        output << "        {" << "\n";
        output << "#<!--================== HisGram ====================-->" << "\n";
        output << "" << "\n";
        output << "            $Sort_Array = Array();" << "\n";
        output << "            $Data_Array = Array();" << "\n";
        output << "" << "\n";
        output << "            $Data_Array[ '5p' ] = Array();" << "\n";
        output << "            $Data_Array[ '3p' ] = Array();" << "\n";
        output << "" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Tail_Array[0] as $Mir )" << "\n";
        output << "                Array_Push( $Sort_Array, $Mir[1] - $Mir[2] );" << "\n";
        output << "" << "\n";
        output << "            Array_Multisort( $Sort_Array, SORT_NUMERIC, SORT_DESC, $Tail_Array[0] );" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Tail_Array[0] as $Mir )" << "\n";
        output << "            {" << "\n";
        output << "                $Data_Array[ '5p' ][ $Mir[0] ] = $Mir[1];" << "\n";
        output << "                $Data_Array[ '3p' ][ $Mir[0] ] = $Mir[2];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>" << "\n";
        output << "                var svg_width  = window.innerWidth;" << "\n";
        output << "                var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "                var margin = {top: 10, right: 20, bottom: 10, left: 20}" << "\n";
        output << "                    width  = svg_width  - margin.left - margin.right," << "\n";
        output << "                    height = svg_height - margin.top  - margin.bottom;" << "\n";
        output << "" << "\n";
        output << "                var svg = d3.select('body').append('div')" << "\n";
        output << "                    .attr('id', 'svg')" << "\n";
        output << "                    .style('width', width + 'px')" << "\n";
        output << "                    .style('height', height + 'px')" << "\n";
        output << "                    .style('display', 'inline-block' )" << "\n";
        output << "                    .append('svg');" << "\n";
        output << "" << "\n";
        output << "                nv.addGraph( function() {" << "\n";
        output << "                    var chart = nv.models.multiBarHorizontalChart()" << "\n";
        output << "                        .x(function(d) { return d.label })" << "\n";
        output << "                        .y(function(d) { return d.value })" << "\n";
        output << "                        .color([ '#d67777', '#4f99b4' ])" << "\n";
        output << "                        .margin({left: 80})" << "\n";
        output << "                        .groupSpacing(0)" << "\n";
        output << "                        .showControls(true)" << "\n";
        output << "                        .showValues(true)" << "\n";
        output << "                        .showXAxis(false)" << "\n";
        output << "                        .stacked(true)" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    d3.select('#svg svg')" << "\n";
        output << "                        .data([data])" << "\n";
        output << "                        .call(chart);" << "\n";
        output << "" << "\n";
        output << "                    nv.utils.windowResize(chart.update);" << "\n";
        output << "                    return chart;" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "                var data = [\";" << "\n";
        output << "" << "\n";
        output << "            Foreach( $Data_Array as $Key => $Data )" << "\n";
        output << "            {" << "\n";
        output << "                echo '{key:\"'.$Key.'\",values:[';" << "\n";
        output << "" << "\n";
        output << "                Foreach( $Data as $Label => $Value )" << "\n";
        output << "                {" << "\n";
        output << "                    echo '{label:\"'.$Label.'\",value:'.$Value.'},';" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                echo ']},';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '];</script>';" << "\n";
        output << "        }" << "\n";
        output << "        else" << "\n";
        output << "        {" << "\n";
        output << "#<!--================== BoxPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "            $Boxs_Array = Array();" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $Tail_Array ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                Sort( $Tail_Array[$i] );" << "\n";
        output << "" << "\n";
        output << "                $Boxs_Array[$i][ 'Q1' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.25 )];" << "\n";
        output << "                $Boxs_Array[$i][ 'Q2' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.5  )];" << "\n";
        output << "                $Boxs_Array[$i][ 'Q3' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * 0.75 )];" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_low'  ] = $Tail_Array[$i][0];" << "\n";
        output << "                $Boxs_Array[$i][ 'whisker_high' ] = $Tail_Array[$i][ Count( $Tail_Array[$i] ) -1];" << "\n";
        output << "" << "\n";
        output << "                if( $isTrimmed != '' && $isTrimmed != 'isTrimmed' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Rate = $isTrimmed / 100;" << "\n";
        output << "                    $Boxs_Array[$i][ 'whisker_low'  ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) *       $Rate  )];" << "\n";
        output << "                    $Boxs_Array[$i][ 'whisker_high' ] = $Tail_Array[$i][ Floor(( Count( $Tail_Array[$i] ) -1 ) * ( 1 - $Rate ))];" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( $Boxs_Array[$i][ 'whisker_low'  ] < $yMin ) $yMin = $Boxs_Array[$i][ 'whisker_low'  ];" << "\n";
        output << "                if( $Boxs_Array[$i][ 'whisker_high' ] > $yMax ) $yMax = $Boxs_Array[$i][ 'whisker_high' ];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $yMin = $DiffType == 'Ratios_5p3p' ? 0 : -1;" << "\n";
        output << "            $yMax = 1;" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>" << "\n";
        output << "                var svg_width  = window.innerWidth;" << "\n";
        output << "                var svg_height = window.innerHeight;" << "\n";
        output << "" << "\n";
        output << "                var margin = {top: 10, right: 20, bottom: 10, left: 20}" << "\n";
        output << "                    width  = svg_width  - margin.left - margin.right," << "\n";
        output << "                    height = svg_height - margin.top  - margin.bottom;" << "\n";
        output << "" << "\n";
        output << "                var svg = d3.select('body').append('div')" << "\n";
        output << "                    .attr('id', 'svg')" << "\n";
        output << "                    .style('width', width + 'px')" << "\n";
        output << "                    .style('height', height + 'px')" << "\n";
        output << "                    .style('display', 'inline-block' )" << "\n";
        output << "                    .append('svg');" << "\n";
        output << "" << "\n";
        output << "                nv.addGraph( function() {" << "\n";
        output << "                    var chart = nv.models.boxPlotChart()" << "\n";
        output << "                        .x(function(d) { return d.label })" << "\n";
        output << "                        .staggerLabels(true)" << "\n";
        output << "                        .maxBoxWidth(75) // prevent boxes from being incredibly wide" << "\n";
        output << "                        .yDomain([$yMin, $yMax]);" << "\n";
        output << "" << "\n";
        output << "                    d3.select('#svg svg')" << "\n";
        output << "                        .data([data])" << "\n";
        output << "                        .call(chart);" << "\n";
        output << "" << "\n";
        output << "                    nv.utils.windowResize(chart.update);" << "\n";
        output << "                    return chart;" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "                var data = [\";" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $Boxs_Array ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '{label:\"'.Explode( '.', $TSVs[$i] )[0].'\",values:{';" << "\n";
        output << "" << "\n";
        output << "                Foreach( $Boxs_Array[$i] as $Q => $Value )" << "\n";
        output << "                    echo $Q.':'.$Value.',';" << "\n";
        output << "" << "\n";
        output << "                echo 'outliers:[]}}';" << "\n";
        output << "                if( $i != Count( $Boxs_Array )-1 ) echo \",\\n\";" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '];</script>';" << "\n";
        output << "        }" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
