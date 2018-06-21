#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerTaildot
{

  public:

    static std::vector< std::map< std::string, std::vector< double >>> anno_tail_table;
    static std::vector< std::map< std::string, std::vector< double >>> isom_tail_table;

    GeneTypeAnalyzerTaildot()
    {}

    static void make_taildot_table(
            const std::string& biotype,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< BedSampleType >& bed_samples,
            auto& genome_table
            )
    {
        anno_tail_table.clear();
        isom_tail_table.clear();

        std::vector< std::string > split;
        std::map< std::string, std::vector< double >> anno_map, isom_map;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            anno_map.emplace( split[0], std::vector< double >( 8, 0.0 ));
            isom_map.emplace( anno    , std::vector< double >( 8, 0.0 ));
            // GMPM TailLens    Tailing％   A_Tail％    C_Tail％    G_Tail％    T_Tail％    Other_Tail％
            // 0    1           2           3           4           5           6           7
        }

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            anno_tail_table.emplace_back( anno_map );
            make_tail_counting_table( biotype, bed_samples[ smp ], anno_tail_table[ smp ], genome_table );

            isom_tail_table.emplace_back( isom_map );
            make_tail_counting_table( biotype, bed_samples[ smp ], isom_tail_table[ smp ], genome_table, true );
        }
    }

    static void make_tail_counting_table(
            const std::string& biotype,
            BedSampleType& bed_sample,
            std::map< std::string, std::vector< double >>& tail_table,
            auto& genome_table,
            bool is_isomir = false
            )
    {
        std::string annotation;
        double length, sum;
        std::size_t tail;

        std::vector< double > anno_tailing_vec;
        std::map< double, double > lens_p;

        std::map< std::string, std::pair< double, double >> anno_gmpm;
        std::map< std::string, std::map< std::size_t, std::map< double, double >>> anno_temp;
        //          annotation              tail                length  ppm

        for( auto& raw_bed : bed_sample.second )
        {
            // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
            {
                std::size_t i = 0; // do first priority
                if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                {
                    if( raw_bed.annotation_info_[i][0] == biotype || ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" )))
                    {
                        tail = GeneTypeAnalyzerCounting::which_tail( raw_bed.getTail() );
                        length = ( double )( raw_bed.length_ );

                        for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                        {
                            annotation = raw_bed.annotation_info_[i][ j+1 ] + ( !is_isomir ? ""
                                : ( "_" + raw_bed.getReadSeq( genome_table ).substr( 1, 7 )
                                + ( raw_bed.seed_md_tag != "" ? ( "|" + raw_bed.seed_md_tag ) : "" )
                                ));


                            if( anno_gmpm.find( annotation ) == anno_gmpm.end() )
                                anno_gmpm[ annotation ] = std::make_pair( 0.0, 0.0 );

                            if( tail == 5 )
                            {
                                anno_gmpm[ annotation ].first += raw_bed.ppm_;
                                continue;
                            }

                            if( anno_temp[ annotation ][ tail ].find( length ) == anno_temp[ annotation ][ tail ].end() )
                                anno_temp[ annotation ][ tail ][ length ] = 0.0;

                            anno_temp[ annotation ][ tail ][ length ] += raw_bed.ppm_;
                            anno_gmpm[ annotation ].second += raw_bed.ppm_;
                        }
                    }
                }
            }
        }

        for( auto& anno : anno_gmpm )
        {
            sum = 0.0;
            length = 0.0;

            lens_p.clear();
            anno_tailing_vec = std::vector< double >( 8, 0.0 );

            for( auto& tail : anno_temp[ anno.first ] ) for( auto& len : tail.second )
            {
                if( lens_p.find( len.first ) == lens_p.end() )
                    lens_p[ len.first ] = 0;

                lens_p[ len.first ] += 1;
                sum += 1;
            }

            for( auto& tail : anno_temp[ anno.first ] ) for( auto& len : tail.second )
            {
                anno_tailing_vec[ tail.first + 3 ] += len.second * ( lens_p[ len.first ] / sum ); 
                len.second = len.second / anno.second.second;
                length += len.first * len.second;
            }

            if( length == 0.0 ) continue;

            anno_tailing_vec[0] = anno.second.first + anno.second.second;   // GMPM
            anno_tailing_vec[1] = length;                                   // TailLens
            anno_tailing_vec[2] = anno.second.second /( anno.second.first + anno.second.second );   // Tailing％
            anno_tailing_vec[3] = anno_tailing_vec[3] / anno.second.second; // A_Tail％
            anno_tailing_vec[4] = anno_tailing_vec[4] / anno.second.second; // C_Tail％
            anno_tailing_vec[5] = anno_tailing_vec[5] / anno.second.second; // G_Tail％
            anno_tailing_vec[6] = anno_tailing_vec[6] / anno.second.second; // T_Tail％
            anno_tailing_vec[7] = anno_tailing_vec[7] / anno.second.second; // Other_Tail％

            tail_table[ anno.first ] = anno_tailing_vec;
        }
    }

    static void output_taildot_isomirs(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name,
            const std::size_t& smp
            )
    {
        std::ofstream output( output_name + sample_name + "-isomiRs.tsv" );
        output << sample_name << "\tGMPM\tTailLens\tTailing％\tA_Tail％\tC_Tail％\tG_Tail％\tT_Tail％\tOther_Tail％\n";

        for( auto& anno : ano_len_idx.first )
        {
            output << anno << ( anno_mark.find( anno ) != anno_mark.end() ? anno_mark[ anno ] : "" );

            for( auto& ratio : isom_tail_table[ smp ][ anno ] )
                output << "\t" << std::setprecision( 2 ) << std::fixed << ratio;

            output << "\n";
        }

        output.close();
    }

    static void output_taildot(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::map< std::string, std::string >& anno_mark,
            const std::string& sample_name,
            const std::size_t& smp
            )
    {
        std::vector< std::string > split;

        std::set< std::string > anno_idx_set;
        std::set< std::string > anno_mark_set;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            anno_idx_set.emplace( split[0] );

            if( anno_mark.find( anno ) != anno_mark.end() )
            {
                if( anno_mark[ anno ] == "!" || anno_mark[ anno ] == "*!" )
                    anno_mark_set.emplace( split[0] );
            }
        }

        std::ofstream output( output_name + sample_name + ".tsv" );
        output << sample_name << "\tGMPM\tTailLens\tTailing％\tA_Tail％\tC_Tail％\tG_Tail％\tT_Tail％\tOther_Tail％\n";

        for( auto& anno : anno_idx_set )
        {
            output << anno << ( anno_mark_set.find( anno ) != anno_mark_set.end() ? "!" : "" );

            for( auto& ratio : anno_tail_table[ smp ][ anno ] )
                output << "\t" << std::setprecision( 2 ) << std::fixed << ratio;

            output << "\n";
        }

        output.close();
    }

    static void output_taildot_visualization( const std::string& output_name )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <body>" << "\n";
        output << "" << "\n";
        output << "    <? " << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "        $RatioType = $_POST['RatioType'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $LenMin = $_POST['LenMin'];" << "\n";
        output << "        $LenMax = $_POST['LenMax'];" << "\n";
        output << "        $IsomiRs = $_POST['IsomiRs'];" << "\n";
        output << "        $RatioMin = $_POST['RatioMin'];" << "\n";
        output << "        $RatioMax = $_POST['RatioMax'];" << "\n";
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        $isAbundant = $_POST['isAbundant'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "" << "\n";
        output << "        echo '<style type=\"text/css\">" << "\n";
        output << "                div[id=\"Dot\"] {" << "\n";
        output << "                    position: absolute; " << "\n";
        output << "                    text-align: left; " << "\n";
        output << "                    width: 120px;" << "\n";
        output << "                    height: 110px; " << "\n";
        output << "                    padding: 8px; " << "\n";
        output << "                    background: lightsteelblue; " << "\n";
        output << "                    font: 12px sans-serif;" << "\n";
        output << "                    border: 0px;" << "\n";
        output << "                    border-radius: 8px; " << "\n";
        output << "                    pointer-events: none; " << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                th {" << "\n";
        output << "                    border: 1px solid black;" << "\n";
        output << "                } " << "\n";
        output << "" << "\n";
        output << "                .axis path, .axis line {" << "\n";
        output << "                    fill: none;" << "\n";
        output << "                    stroke: #000;" << "\n";
        output << "                    shape-rendering: crispEdges;" << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                #grid {" << "\n";
        output << "                    position: fixed;" << "\n";
        output << "                    width: 100%;" << "\n";
        output << "                    bottom: 0;" << "\n";
        output << "                    height: 300px;" << "\n";
        output << "                }" << "\n";
        output << "                " << "\n";
        output << "                .slick-row:hover {" << "\n";
        output << "                    font-weight: bold;" << "\n";
        output << "                    color: #069;" << "\n";
        output << "                }" << "\n";
        output << "            </style>';" << "\n";
        output << "" << "\n";
        output << "#<!--================ RatioType =================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=RatioType onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($RatioType=='') echo 'selected'; echo '>RatioType</option>';" << "\n";
        output << "" << "\n";
        output << "        $RatioType_List = array('Tailing％', 'A_Tail％', 'C_Tail％', 'G_Tail％', 'T_Tail％', 'Other_Tail％');" << "\n";
        output << "" << "\n";
        output << "        $RatioType_Size = Count( $RatioType_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $RatioType_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$RatioType_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $RatioType == $RatioType_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $RatioType_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='LenMin' value='$LenMin' />" << "\n";
        output << "            <input type='hidden' name='LenMax' value='$LenMax' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='RatioMin' value='$RatioMin' />" << "\n";
        output << "            <input type='hidden' name='RatioMax' value='$RatioMax' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
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
        output << "            <input type='hidden' name='RatioType' value='$RatioType' />" << "\n";
        output << "            <input type='hidden' name='LenMin' value='$LenMin' />" << "\n";
        output << "            <input type='hidden' name='LenMax' value='$LenMax' />" << "\n";
        output << "            <input type='hidden' name='RatioMin' value='$RatioMin' />" << "\n";
        output << "            <input type='hidden' name='RatioMax' value='$RatioMax' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "        $TSV_List_Temp = array();" << "\n";
        output << "" << "\n";
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
        output << "" << "\n";
        output << "        $TSV_List = $TSV_List_Temp;" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
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
        output << "        echo '</select>';" << "\n";
        output << "" << "\n";
        output << "        $TSV_File_Temp = $TSV_File.( $IsomiRs == 'Yes' ? '-isomiRs.tsv' : '.tsv' );" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='RatioType' value='$RatioType' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='LenMin' value='$LenMin' />" << "\n";
        output << "            <input type='hidden' name='LenMax' value='$LenMax' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='RatioMin' value='$RatioMin' />" << "\n";
        output << "            <input type='hidden' name='RatioMax' value='$RatioMax' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== is_Abundant ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n";
        output << "" << "\n";
        output << "        $isAbundant_List = array('MostAbundant', 'AllmiRNA');" << "\n";
        output << "        $isAbundant_Size = Count( $isAbundant_List );" << "\n";
        output << "" << "\n";
        output << "        if( $isAbundant == '' )" << "\n";
        output << "            $isAbundant = 'MostAbundant';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$isAbundant_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $isAbundant == $isAbundant_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $isAbundant_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='RatioType' value='$RatioType' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='LenMin' value='$LenMin' />" << "\n";
        output << "            <input type='hidden' name='LenMax' value='$LenMax' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='RatioMin' value='$RatioMin' />" << "\n";
        output << "            <input type='hidden' name='RatioMax' value='$RatioMax' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== ForceMin&Max ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=LenMin size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $LenMin=='' )" << "\n";
        output << "            echo 'LenMin';" << "\n";
        output << "        else" << "\n";
        output << "            echo $LenMin;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        echo '<input type=text name=LenMax size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $LenMax=='' )" << "\n";
        output << "            echo 'LenMax';" << "\n";
        output << "        else" << "\n";
        output << "            echo $LenMax;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        echo '<input type=text name=RatioMin size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $RatioMin=='' )" << "\n";
        output << "            echo 'RatioMin';" << "\n";
        output << "        else" << "\n";
        output << "            echo $RatioMin;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "        echo '<input type=text name=RatioMax size=5 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $RatioMax=='' )" << "\n";
        output << "            echo 'RatioMax';" << "\n";
        output << "        else" << "\n";
        output << "            echo $RatioMax;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter NU ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Filter size=7 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Filter=='' )" << "\n";
        output << "            echo 'FilterGMPM';" << "\n";
        output << "        else" << "\n";
        output << "            echo $Filter;" << "\n";
        output << "" << "\n";
        output << "        echo \" onfocus=\\\"{this.value='';}\\\">" << "\n";
        output << "            <input type='hidden' name='RatioType' value='$RatioType' />" << "\n";
        output << "            <input type='hidden' name='IsomiRs' value='$IsomiRs' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form><br/>\";" << "\n";
        output << "" << "\n";
        output << "#<!--===================== Read File ======================-->" << "\n";
        output << "" << "\n";
        output << "        $Header = Array();" << "\n";
        output << "" << "\n";
        output << "        $AColumn = 0;" << "\n";
        output << "        $RColumn = 0;" << "\n";
        output << "        $TColumn = 1;" << "\n";
        output << "        $LColumn = 2;" << "\n";
        output << "" << "\n";
        output << "        $minLenth = 30;" << "\n";
        output << "        $maxLenth = 0;" << "\n";
        output << "        $minRatio = 0;" << "\n";
        output << "        $maxRatio = 0;" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( $TSV_File_Temp );" << "\n";
        output << "        $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "        $Temp = Tempnam( '/tmp', $RatioType.'_'.$TSV_File_Temp.'_'.$isAbundant.'_'.$Filter );" << "\n";
        output << "        $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n";
        output << "            $miRNA_Seed  = Explode( '*', $inFile_Line[ $AColumn ]);" << "\n";
        output << "" << "\n";
        output << "            if( $i == 0 )" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $inFile_Line ); ++$j )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $inFile_Line[$j] == $RatioType )" << "\n";
        output << "                        $RColumn = $j;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                Array_Push( $Header, Substr( $inFile_Line[ $AColumn ], 0, Strlen( $inFile_Line[ $AColumn ]) ));" << "\n";
        output << "                Array_Push( $Header, Substr( $inFile_Line[ $LColumn ], 0, Strlen( $inFile_Line[ $LColumn ]) ));" << "\n";
        output << "                Array_Push( $Header, Substr( $inFile_Line[ $RColumn ], 0, Strlen( $inFile_Line[ $RColumn ]) -3 ));" << "\n";
        output << "                Array_Push( $Header, Substr( $inFile_Line[ $TColumn ], 0, Strlen( $inFile_Line[ $TColumn ]) ));" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"miRNA\\t\".$Header[1].\"\\t\".$Header[2].\"\\t\".$Header[3].\"\\n\" );" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "            {" << "\n";
        output << "                if( $isAbundant == 'MostAbundant' && $IsomiRs == 'Yes' && Count( $miRNA_Seed ) == 1 ) continue;" << "\n";
        output << "                if( $Filter != 'FilterGMPM' && $inFile_Line[ $TColumn ] < $Filter ) continue;" << "\n";
        output << "" << "\n";
        output << "                if( $minLenth > $inFile_Line[ $LColumn ] && $inFile_Line[ $LColumn ] != 0 ) $minLenth = $inFile_Line[ $LColumn ];" << "\n";
        output << "                if( $maxLenth < $inFile_Line[ $LColumn ] ) $maxLenth = $inFile_Line[ $LColumn ];" << "\n";
        output << "                if( $maxRatio < $inFile_Line[ $RColumn ] ) $maxRatio = $inFile_Line[ $RColumn ];" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp," << "\n";
        output << "                    $miRNA_Seed[0].\"\\t\"." << "\n";
        output << "                    $inFile_Line[ $LColumn ].\"\\t\"." << "\n";
        output << "                    $inFile_Line[ $RColumn ].\"\\t\"." << "\n";
        output << "                    $inFile_Line[ $TColumn ].\"\\n\"" << "\n";
        output << "                );" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "        $minLenth = ( $LenMin   != 'LenMin'   && $LenMin   != '' ? $LenMin   : $minLenth );" << "\n";
        output << "        $minRatio = ( $RatioMin != 'RatioMin' && $RatioMin != '' ? $RatioMin : $minRatio );" << "\n";
        output << "" << "\n";
        output << "        $maxLenth = ( $LenMax   != 'LenMax'   && $LenMax   != '' ? $LenMax   : $maxLenth );" << "\n";
        output << "        $maxRatio = ( $RatioMax != 'RatioMax' && $RatioMax != '' ? $RatioMax : $maxRatio );" << "\n";
        output << "" << "\n";
        output << "#<!--================== DotPlot ====================-->" << "\n";
        output << "" << "\n";
        output << "        $miR_End = ( $IsomiRs == 'Yes' ? 11 : 3 );" << "\n";
        output << "" << "\n";
        output << "        echo \"<script>" << "\n";
        output << "            var svg_width  = window.innerWidth;" << "\n";
        output << "            var svg_height = window.innerHeight;" << "\n";
        output << "        " << "\n";
        output << "            var margin = {top: 20, right: 60, bottom: 70, left: 60}" << "\n";
        output << "                width = svg_width - margin.left - margin.right," << "\n";
        output << "                height = svg_height - margin.top - margin.bottom;" << "\n";
        output << "        " << "\n";
        output << "            var x = d3.scale.linear()" << "\n";
        output << "                .range([0, width]);" << "\n";
        output << "        " << "\n";
        output << "            var y = d3.scale.linear()" << "\n";
        output << "                .range([height, 0]);" << "\n";
        output << "        " << "\n";
        output << "            var color = d3.scale.category10();" << "\n";
        output << "        " << "\n";
        output << "            var xAxis = d3.svg.axis()" << "\n";
        output << "                .scale(x)" << "\n";
        output << "                .orient('bottom');" << "\n";
        output << "        " << "\n";
        output << "            var yAxis = d3.svg.axis()" << "\n";
        output << "                .scale(y)" << "\n";
        output << "                .orient('left');" << "\n";
        output << "        " << "\n";
        output << "            var svg = d3.select('body').append('svg')" << "\n";
        output << "                .attr('width', width+70)" << "\n";
        output << "                .attr('height', height+60)" << "\n";
        output << "            .append('g')" << "\n";
        output << "                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n";
        output << "        " << "\n";
        output << "            d3.tsv('$Temp', function(error, data) {" << "\n";
        output << "        " << "\n";
        output << "                var xExtent = d3.extent(data, function(d) { return d.$Header[2]; });" << "\n";
        output << "                var yExtent = d3.extent(data, function(d) { return d.$Header[1]; });" << "\n";
        output << "        " << "\n";
        output << "                xExtent[0] = $minLenth;" << "\n";
        output << "                yExtent[0] = $minRatio;" << "\n";
        output << "        " << "\n";
        output << "                xExtent[1] = $maxLenth;" << "\n";
        output << "                yExtent[1] = $maxRatio;" << "\n";
        output << "        " << "\n";
        output << "                x.domain( xExtent ).nice();" << "\n";
        output << "                y.domain( yExtent ).nice();" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'x axis')" << "\n";
        output << "                    .attr('transform', 'translate(0,' + height + ')')" << "\n";
        output << "                    .call(xAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('x', width)" << "\n";
        output << "                    .attr('y', -6)" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('$Header[1]');" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'y axis')" << "\n";
        output << "                    .call(yAxis)" << "\n";
        output << "                    .append('text')" << "\n";
        output << "                    .attr('class', 'label')" << "\n";
        output << "                    .attr('transform', 'rotate(-90)')" << "\n";
        output << "                    .attr('y', 6)" << "\n";
        output << "                    .attr('dy', '.71em')" << "\n";
        output << "                    .style('text-anchor', 'end')" << "\n";
        output << "                    .text('$Header[2]');" << "\n";
        output << "        " << "\n";
        output << "                svg.append('g')" << "\n";
        output << "                    .attr('class', 'xyline')" << "\n";
        output << "                    .append('line')" << "\n";
        output << "                    .attr('x1', x($minLenth) )" << "\n";
        output << "                    .attr('y1', y($minRatio) )" << "\n";
        output << "                    .attr('x2', x($maxLenth) )" << "\n";
        output << "                    .attr('y2', y($maxRatio) )" << "\n";
        output << "                    .style('stroke','rgb(255,0,0)')" << "\n";
        output << "                    .style('stroke-width','1');" << "\n";
        output << "        " << "\n";
        output << "                var div = d3.select('body').append('div')" << "\n";
        output << "                    .attr('class', 'tooltip')" << "\n";
        output << "                    .attr('id', 'Dot')" << "\n";
        output << "                    .style('opacity', 0);" << "\n";
        output << "        " << "\n";
        output << "                svg.selectAll('.dot')" << "\n";
        output << "                    .data(data)" << "\n";
        output << "                    .enter()" << "\n";
        output << "                    // .append('a')" << "\n";
        output << "                    // .attr('xlink:href', function(d){ return '../SqAlign/index.php?TSV_File=$TSV_File1.tsv&Annotation_Select=' + d.miRNA.substring( 0, d.miRNA.length - $miR_End )})" << "\n";
        output << "                    // .attr('target', '_blank')" << "\n";
        output << "                    .append('circle')" << "\n";
        output << "                    .attr('class', 'dot')" << "\n";
        output << "                    .attr('r', 5)" << "\n";
        output << "                    .attr('cx', function(d) { return x(d.$Header[1]); })" << "\n";
        output << "                    .attr('cy', function(d) { return y(d.$Header[2]); })" << "\n";
        output << "                    .style('fill', function(d) { return color(d.species); })" << "\n";
        output << "                    .on('mouseover', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(200)" << "\n";
        output << "                            .style('opacity', .9);" << "\n";
        output << "        " << "\n";
        output << "                        div.html('<table align=center ><tr><th colspan=2 >' + d.miRNA +" << "\n";
        output << "                                 '</th></tr><tr><th>Sample</th><th>$Header[0]</th></tr><tr><th>TotalPPM</th><th>' +" << "\n";
        output << "                                 d.$Header[3] + '</th></tr><tr><th>$Header[1]</th><th>' +" << "\n";
        output << "                                 d.$Header[1] + '</th></tr><tr><th>$Header[2]</th><th>' +" << "\n";
        output << "                                 d.$Header[2] + '</th></tr></table>')" << "\n";
        output << "                            .style('left', (d3.event.pageX) + 'px')" << "\n";
        output << "                            .style('top', (d3.event.pageY - 28) + 'px');" << "\n";
        output << "                    })" << "\n";
        output << "                    .on('mouseout', function(d) {" << "\n";
        output << "                        div.transition()" << "\n";
        output << "                            .duration(500)" << "\n";
        output << "                            .style('opacity', 0);" << "\n";
        output << "                    });" << "\n";
        output << "        " << "\n";
        output << "                    window.onresize = function(){" << "\n";
        output << "                    };" << "\n";
        output << "            });" << "\n";
        output << "        </script>\";" << "\n";
        output << "        " << "\n";
        output << "        " << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

std::vector< std::map< std::string, std::vector< double >>> GeneTypeAnalyzerTaildot::anno_tail_table;
std::vector< std::map< std::string, std::vector< double >>> GeneTypeAnalyzerTaildot::isom_tail_table;

} // end of namespace algorithm
} // end of namespace ago
