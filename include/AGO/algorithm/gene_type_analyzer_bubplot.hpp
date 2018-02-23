#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBubplot
{
    using ChrRangeType = std::tuple< std::string, std::size_t, std::size_t, char >;

  public:

    GeneTypeAnalyzerBubplot()
    {}

    static void output_bubplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number,
            auto& genome_table
            )
    {
        for( auto& smp : bed_samples ) if( !boost::filesystem::exists( output_name + smp.first + ".tsv" ))
            boost::filesystem::create_symlink(( "../LenDist/" + smp.first + ".tsv" ).c_str(), ( output_name + smp.first + ".tsv" ).c_str() );

        std::map< std::string, ChrRangeType > chr_mapping = get_chrmap_table( bed_samples, biotype, thread_number );
        std::vector< std::string > out_vec = sequence_formating( chr_mapping, genome_table );

        std::ofstream output( output_name + "AnnoSeq.tsv" );
        output << "Annotation\t5P\t3P";

        for( auto& res : out_vec )
            output << res;

        output << "\n";
        output.close();
    }

    static std::map< std::string, ChrRangeType > get_chrmap_table(
            const std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number
            )
    {
        std::map< std::string, std::vector< std::pair< ChrRangeType, std::size_t >>> chr_mappings;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                for( auto& raw_bed_info : raw_bed.annotation_info_ )
                {
                    for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                    {
                        if( raw_bed_info[i] != biotype ) continue;

                        if( chr_mappings.find( raw_bed_info[ i+1 ] ) == chr_mappings.end() )
                            chr_mappings[ raw_bed_info[ i+1 ]] = std::vector< std::pair< ChrRangeType, std::size_t >>();

                        chr_mappings[ raw_bed_info[ i+1 ]].emplace_back(
                                std::make_pair( std::make_tuple(
                                    raw_bed.chromosome_,
                                    raw_bed.start_,
                                    raw_bed.end_,
                                    raw_bed.strand_ ),
                                    0
                                ));
                    }
                }
            }
        }

        ParaThreadPool parallel_pool( thread_number );
        std::vector< std::pair< std::string, std::vector< std::pair< ChrRangeType, std::size_t >>>> parallel_vec;
        std::vector< std::size_t > parallel_indx;

        std::size_t task_number = thread_number * 10;

        for( auto& anno : chr_mappings )
            parallel_vec.emplace_back( anno );

        chr_mappings.clear();

        for( std::size_t anno = 0; anno < parallel_vec.size(); ++anno )
        {
            parallel_indx.emplace_back( anno );

            if( parallel_indx.size() >= task_number )
                parallel_pool.job_post([ parallel_indx, &parallel_vec ] ()
                {
                    for( auto& idx : parallel_indx )
                    {
                        recursive_merge( parallel_vec[ idx ].second, 0 );
                        range_counting_sort( parallel_vec[ idx ].second );
                    }
                });

            parallel_indx.clear();
        }

        parallel_pool.flush_pool();

        std::map< std::string, ChrRangeType > chr_mapping_res;

        for( auto& anno : parallel_vec )
            chr_mapping_res[ anno.first ] = anno.second[0].first;

        return chr_mapping_res;
    }

    static void recursive_merge( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges, std::size_t start_idx )
    {
        if( start_idx == ranges.size() ) return;
        std::size_t counts  = ranges[ start_idx ].second;

        for( std::size_t idx = 0; idx < ranges.size(); ++idx )
        {
            if( idx == start_idx ) continue;
            
            if( std::get<0>( ranges[ idx ].first ) != std::get<0>( ranges[ start_idx ].first )) continue;
            if( std::get<3>( ranges[ idx ].first ) != std::get<3>( ranges[ start_idx ].first )) continue;

            if( std::get<1>( ranges[ idx ].first ) <= std::get<2>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) >= std::get<1>( ranges[ start_idx ].first ) || 
                std::get<1>( ranges[ idx ].first ) >= std::get<1>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) <= std::get<2>( ranges[ start_idx ].first ) || 
                std::get<1>( ranges[ idx ].first ) <= std::get<1>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) >= std::get<2>( ranges[ start_idx ].first ) )
            {
                if( std::get<1>( ranges[ idx ].first ) < std::get<1>( ranges[ start_idx ].first )) std::get<1>( ranges[ start_idx ].first ) = std::get<1>( ranges[ idx ].first );
                if( std::get<2>( ranges[ idx ].first ) > std::get<2>( ranges[ start_idx ].first )) std::get<2>( ranges[ start_idx ].first ) = std::get<2>( ranges[ idx ].first );
                if( idx < start_idx ) start_idx--;

                ranges.erase( ranges.begin(), ranges.begin() + idx + 1 );
                counts++;
                idx--;
            }
        }

        if( counts == ranges[ start_idx ].second ) recursive_merge( ranges, start_idx +1 );
    }

    static void range_counting_sort( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges )
    {
        std::sort( ranges.begin(), ranges.end(),
            []( const std::pair< ChrRangeType, std::size_t >& a, const std::pair< ChrRangeType, std::size_t >& b )
            { return a.second > b.second; });
    }

    static std::vector< std::string > sequence_formating(
            const std::map< std::string, ChrRangeType >& chr_mapping,
            auto& genome_table
            )
    {
        std::vector< std::string > res_vec;
        std::map< std::string, std::tuple< std::string, std::string >> temp_map;

        for( auto& anno : chr_mapping )
        {
            if( temp_map.find( anno.first.substr( 0, anno.first.length() -3 )) == temp_map.end() )
                temp_map[ anno.first.substr( 0, anno.first.length() -3 )] = std::tuple< std::string, std::string >();

            switch( anno.first.at( anno.first.length() -2 ))
            {
                case '5' : std::get<0>( temp_map[ anno.first.substr( 0, anno.first.length() -3 )]) = get_sequence( anno.second, genome_table ); break;
                case '3' : std::get<1>( temp_map[ anno.first.substr( 0, anno.first.length() -3 )]) = get_sequence( anno.second, genome_table ); break;
            }
        }

        for( auto& anno : temp_map )
        {
            res_vec.emplace_back(
                    "\n" + anno.first +
                    "\t" + std::get<0>( anno.second ) +
                    "\t" + std::get<1>( anno.second ) );
        }

        return res_vec;
    }

    static std::string get_sequence( const ChrRangeType& range, auto& genome_table )
    {
		std::string read_seq =
            genome_table[ std::get<0>( range )].substr( std::get<1>( range ) - 1, std::get<2>( range ) - std::get<1>( range ));

		std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper );

		if( std::get<3>( range ) == '-' )
		{
			std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), []( char c ){ return complement( c ); });
			std::reverse( read_seq.begin(), read_seq.end() );
		}

        return read_seq;
    }

	static char complement( char c )
	{
		switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
		}
		return c;
	}

    static void output_bubplot_visualization( const std::string& output_name )
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
        output << "        $Chart_Types = $_POST['Chart_Types'];" << "\n";
        output << "        $GMPMT_Types = $_POST['GMPMT_Types'];" << "\n";
        output << "" << "\n";
        output << "        $isLog2 = $_POST['isLog2'];" << "\n";
        output << "        $Annotation_Select = $_POST['Annotation_Select'];" << "\n";
        output << "" << "\n";
        output << "        $Annotation_Arms = $_POST['Annotation_Arms'];" << "\n";
        output << "        $Annotation_Arms_Array = $_POST['Annotation_Arms_Array'];" << "\n";
        output << "        $isArmSelect = $_POST['isArmSelect'];" << "\n";
        output << "" << "\n";
        output << "        $Sample_Files = $_POST['Sample_Files'];" << "\n";
        output << "        $Sample_Files_Array = $_POST['Sample_Files_Array'];" << "\n";
        output << "        $isSampleSelect = $_POST['isSampleSelect'];" << "\n";
        output << "" << "\n";
        output << "#<!--================== ChartTypes Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Chart_Types onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Chart_Types=='') echo 'selected'; echo '>Charts</option>';" << "\n";
        output << "" << "\n";
        output << "        $Charts_List = array('bubble', 'heatmap');" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Charts_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Charts_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Chart_Types == $Charts_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $Charts_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Arms' value='$Annotation_Arms' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPMT_Types Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=GMPMT_Types onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($GMPMT_Types=='') echo 'selected'; echo '>Types</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPMT_List = array('GMPM', 'GM', 'PM', 'A_Tail', 'C_Tail', 'G_Tail', 'T_Tail', 'Other_Tail');" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $GMPMT_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$GMPMT_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $GMPMT_Types == $GMPMT_List[$i] )" << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>' . $GMPMT_List[$i] . '</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Arms' value='$Annotation_Arms' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( './AnnoSeq.tsv' );" << "\n";
        output << "        $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "        $Arm_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "        {" << "\n";
        output << "            $anno = Explode( \"\t\", $inFile_Lines[$j] );" << "\n";
        output << "" << "\n";
        output << "            Array_Push( $Anno_Array, $anno[0] );" << "\n";
        output << "" << "\n";
        output << "            if( $Annotation_Select == $anno[0] )" << "\n";
        output << "            {" << "\n";
        output << "                for( $k = 1; $k < Count( $anno ); ++$k )" << "\n";
        output << "                {" << "\n";
        output << "                    Array_Push( $Arm_Array, $anno[$k] );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Sort( $Anno_Array, SORT_STRING );" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Annotation_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Anno_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Anno_Array[$i] != '' )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$Anno_Array[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $Annotation_Select == $Anno_Array[$i] ) " << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "                echo '>'.$Anno_Array[$i].'</option>';" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Arms' value='$Annotation_Arms' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== AnnotatArm Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $Arms_List = array('5p', '3p');" << "\n";
        output << "        if( $isArmSelect ) $Annotation_Arms = Implode( \",\", $Annotation_Arms_Array );" << "\n";
        output << "        $Annotation_Arms_Array = Explode( \",\", $Annotation_Arms );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Arm_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=checkbox onchange=this.form.submit(); name=Annotation_Arms_Array[] value='.$Arms_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            For( $j = 0; $j < Count( $Annotation_Arms_Array ); ++$j )" << "\n";
        output << "                if( $Annotation_Arms_Array[$j] == $Arms_List[$i] )" << "\n";
        output << "                    echo 'checked ';" << "\n";
        output << "" << "\n";
        output << "            if( $Arm_Array[$i] == '' )" << "\n";
        output << "                echo 'disabled=\"disabled\" ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Arms_List[$i].'</input>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='isArmSelect' value='true' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== SampleFile Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv | grep -v AnnoSeq' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        if( $isSampleSelect ) $Sample_Files = Implode( \",\", $Sample_Files_Array );" << "\n";
        output << "        $Sample_Files_Array = Explode( \",\", $Sample_Files );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=checkbox onchange=this.form.submit(); name=Sample_Files_Array[] value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            For( $j = 0; $j < Count( $Sample_Files_Array ); ++$j )" << "\n";
        output << "                if( $Sample_Files_Array[$j] == $TSV_List[$i] )" << "\n";
        output << "                    echo 'checked ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.Explode( \".\", $TSV_List[$i] )[0].'</input>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Arms' value='$Annotation_Arms' />" << "\n";
        output << "            <input type='hidden' name='isSampleSelect' value='true' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--======================== isLog2 =========================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=checkbox onchange=this.form.submit(); name=isLog2 value=\"--log2\" ';" << "\n";
        output << "" << "\n";
        output << "        if( $isLog2 == \"--log2\" ) echo 'checked ';" << "\n";
        output << "        echo '>'.Log2.'</input>';" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Arms' value='$Annotation_Arms' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--======================== Debug ==========================-->" << "\n";
        output << "" << "\n";
        output << "        echo \"<br/> node heatmap_bubble_plot.js\";" << "\n";
        output << "        if( $Chart_Types != '' ) echo \" --mode $Chart_Types\";" << "\n";
        output << "        if( $GMPMT_Types != '' ) echo \" --type $GMPMT_Types\";" << "\n";
        output << "        if( $Annotation_Select != '' ) echo \" --input $Annotation_Select\";" << "\n";
        output << "        if( $Annotation_Arms != '' ) echo \" --arm \".Str_Replace( ',', '', $Annotation_Arms );" << "\n";
        output << "        if( $Sample_Files != '' ) echo \"  --files \".Str_Replace( ',', ' ', $Sample_Files );" << "\n";
        output << "        if( $isLog2 != '' ) echo \" $isLog2\";" << "\n";
        output << "" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
