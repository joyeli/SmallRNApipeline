#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBubplot
{
    using ChrRangeType = std::tuple< std::string, std::size_t, std::size_t, char, std::size_t >;
    //                                  chr         start           end   strand   seed_end

  public:

    GeneTypeAnalyzerBubplot()
    {}

    static void output_bubplot(
            const std::string& output_name,
            std::vector< BedSampleType >& bed_samples,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            const AnnoLengthIndexType& ano_len_idx,
            const std::string& biotype,
            const std::size_t& thread_number,
            const std::size_t& extend_merge,
            const std::size_t& filter_ppm,
            auto& genome_table
            )
    {
        std::map< std::string, double > ppm_table = get_ppm_table( ano_len_idx, anno_table_tail );
        std::map< std::string, ChrRangeType > chr_mapping = get_chrmap_table( bed_samples, biotype, thread_number, extend_merge, filter_ppm );
        std::vector< std::string > out_vec = sequence_formating( chr_mapping, ppm_table, genome_table );

        std::ofstream output( output_name + "AnnoSeq.tsv" );
        output << "Annotation\t5P\t3P";

        for( auto& res : out_vec )
            output << res;

        output << "\n";
        output.close();
    }
    
    static std::map< std::string, double > get_ppm_table(
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
    {
        std::vector< std::string > split; 
        std::map< std::string, double > ppm_table;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            if( split.size() > 2 ) for( std::size_t i = 1; i < split.size() -2; ++i )
                split[0] = split[0] + "_" + split[i];

            split[0] = split[0].substr( 0, split[0].length() -3 );

            if( ppm_table.find( split[0] ) == ppm_table.end() )
                ppm_table[ split[0] ] = 0.0;

            for( auto& anno_table : anno_table_tail )
            {
                for( std::size_t i = 0; i < 6; i++ )
                {
                    if( anno_table[i].find( anno ) == anno_table[i].end() )
                        continue;

                    for( auto& len : ano_len_idx.second )
                        if( anno_table[i][ anno ].find( len ) != anno_table[i][ anno ].end() )
                            ppm_table[ split[0] ] += anno_table[i][ anno ][ len ];
                }
            }
        }

        for( auto& anno : ppm_table )
            anno.second = anno.second / (double)( anno_table_tail.size() );

        return ppm_table;
    }

    static std::map< std::string, ChrRangeType > get_chrmap_table(
            std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number,
            const std::size_t& extend_merge,
            const std::size_t& filter_ppm
            )
    {
        ChrRangeType range_temp;
        std::map< std::string, std::map< ChrRangeType, std::size_t >> chr_mappings;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
                {
                    std::size_t i = 0; // do first priority
                    if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                    {
                        if(( raw_bed.annotation_info_[i][0] == biotype ) ||
                           ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))) 
                        {
                            for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                            {
                                range_temp = std::make_tuple(
                                        raw_bed.chromosome_,
                                        ( raw_bed.strand_ == '+' ? raw_bed.start_ : raw_bed.start_ - extend_merge ),
                                        ( raw_bed.strand_ == '+' ? raw_bed.end_ + extend_merge : raw_bed.end_ ),
                                        raw_bed.strand_,
                                        ( raw_bed.strand_ == '+' ? raw_bed.start_ + 8 : raw_bed.end_ - 8 )
                                        );

                                if( chr_mappings[ raw_bed.annotation_info_[i][ j+1 ]].find( range_temp ) == chr_mappings[ raw_bed.annotation_info_[i][ j+1 ]].end() )
                                    chr_mappings[ raw_bed.annotation_info_[i][ j+1 ]][ range_temp ] = 0;

                                chr_mappings[ raw_bed.annotation_info_[i][ j+1 ]][ range_temp ]++;
                            }
                        }
                    }
                }
            }
        }

        ParaThreadPool parallel_pool( thread_number );
        std::size_t task_number = thread_number * 10;

        std::vector< std::size_t > parallel_indx;
        std::vector< std::pair< ChrRangeType, std::size_t >> ranges_temp;
        std::vector< std::pair< std::string, std::vector< std::pair< ChrRangeType, std::size_t >>>> parallel_vec;

        for( auto& anno : chr_mappings )
        {
            for( auto& ranges : anno.second ) ranges_temp.emplace_back( ranges );
            parallel_vec.emplace_back( std::make_pair( anno.first, ranges_temp ));
            ranges_temp.clear();
        }

        chr_mappings.clear();

        for( std::size_t anno = 0; anno < parallel_vec.size(); ++anno )
        {
            parallel_indx.emplace_back( anno );

            if( parallel_indx.size() >= task_number )
            {
                parallel_pool.job_post([ parallel_indx, &parallel_vec ] ()
                {
                    for( auto& idx : parallel_indx )
                    {
                        bool isbreak = false;
                        recursive_merge( parallel_vec[ idx ].second, 0, isbreak );
                        range_counting_sort( parallel_vec[ idx ].second );
                    }
                });

                parallel_indx.clear();
            }
        }

        parallel_pool.job_post([ parallel_indx, &parallel_vec ] ()
        {
            for( auto& idx : parallel_indx )
            {
                bool isbreak = false;
                recursive_merge( parallel_vec[ idx ].second, 0, isbreak );
                range_counting_sort( parallel_vec[ idx ].second );
            }
        });

        parallel_indx.clear();
        parallel_pool.flush_pool();

        std::map< std::string, ChrRangeType > chr_mapping_res;

        for( auto& anno : parallel_vec )
            chr_mapping_res[ anno.first ] = anno.second[0].first;

        return chr_mapping_res;
    }

    static void recursive_merge( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges, std::size_t start_idx, bool& isbreak )
    {
        if( start_idx+1 == ranges.size() ) 
        {
            isbreak = true;
            return;
        }

        std::size_t counts = ranges[ start_idx ].second;

        for( std::size_t idx = 0; idx < ranges.size(); ++idx )
        {
            if( idx == start_idx ) continue;

            if( std::get<0>( ranges[ start_idx ].first ) != std::get<0>( ranges[ idx ].first )) continue;
            if( std::get<3>( ranges[ start_idx ].first ) != std::get<3>( ranges[ idx ].first )) continue;
            
            if(( std::get<1>( ranges[ start_idx ].first ) >= std::get<1>( ranges[ idx ].first ) && std::get<1>( ranges[ start_idx ].first ) <= std::get<2>( ranges[ idx ].first )) ||
               ( std::get<2>( ranges[ start_idx ].first ) <= std::get<2>( ranges[ idx ].first ) && std::get<2>( ranges[ start_idx ].first ) >= std::get<1>( ranges[ idx ].first )) ||
               ( std::get<1>( ranges[ start_idx ].first ) <= std::get<1>( ranges[ idx ].first ) && std::get<2>( ranges[ start_idx ].first ) >= std::get<2>( ranges[ idx ].first )) ||
               ( std::get<1>( ranges[ start_idx ].first ) >= std::get<1>( ranges[ idx ].first ) && std::get<2>( ranges[ start_idx ].first ) <= std::get<2>( ranges[ idx ].first )) ) 
            {
                if( std::get<1>( ranges[ start_idx ].first ) > std::get<1>( ranges[ idx ].first ))
                    std::get<1>( ranges[ start_idx ].first ) = std::get<1>( ranges[ idx ].first );

                if( std::get<2>( ranges[ start_idx ].first ) < std::get<2>( ranges[ idx ].first ))
                    std::get<2>( ranges[ start_idx ].first ) = std::get<2>( ranges[ idx ].first );

                switch( std::get<3>( ranges[ start_idx ].first ))
                {
                    case '+' :
                    if( std::get<4>( ranges[ start_idx ].first ) < std::get<4>( ranges[ idx ].first ))
                        std::get<4>( ranges[ start_idx ].first ) = std::get<4>( ranges[ idx ].first );
                    break;

                    case '-' :
                    if( std::get<4>( ranges[ start_idx ].first ) > std::get<4>( ranges[ idx ].first ))
                        std::get<4>( ranges[ start_idx ].first ) = std::get<4>( ranges[ idx ].first );
                    break;
                }

                ranges[ start_idx ].second += ranges[ idx ].second;
                if( idx < start_idx ) start_idx--;

                ranges.erase( ranges.begin() + idx );
                idx--;
            }
        }

        if( counts == ranges[ start_idx ].second ) recursive_merge( ranges, start_idx +1, isbreak );
        if( !isbreak ) recursive_merge( ranges, start_idx, isbreak );
    }

    static void range_counting_sort( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges )
    {
        std::sort( ranges.begin(), ranges.end(),
            []( const std::pair< ChrRangeType, std::size_t >& a, const std::pair< ChrRangeType, std::size_t >& b )
            { return a.second > b.second; });
    }

    static std::vector< std::string > sequence_formating(
            const std::map< std::string, ChrRangeType >& chr_mapping,
            std::map< std::string, double >& ppm_table,
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
                    "\t" + std::to_string( ppm_table[ anno.first ]) +
                    "\t" + std::get<0>( anno.second ) +
                    "\t" + std::get<1>( anno.second ) );
        }

        return res_vec;
    }

    static std::string get_sequence( const ChrRangeType& range, auto& genome_table )
    {
		std::string read_seq;

        switch( std::get<3>( range ))
        {
            case '+' : read_seq = GeneTypeAnalyzerCounting::seqT2U( genome_table[ std::get<0>( range )].substr( std::get<1>( range )   , std::get<4>( range ) - std::get<1>( range ) -1 )); break;
            case '-' : read_seq = GeneTypeAnalyzerCounting::seqT2U( genome_table[ std::get<0>( range )].substr( std::get<4>( range ) -1, std::get<2>( range ) - std::get<4>( range ) -1 )); break;
        }

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
			case 'A': c = 'U'; break;
			case 'U': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
		}
		return c;
	}

    static void make_file_link( const std::string& p )
    {
        std::set< std::string > samples;
        std::vector< std::string > temp;
        boost::filesystem::path path( p + "../LenDist/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            boost::iter_split( temp, file.path().filename().string(), boost::algorithm::first_finder( "." ));
            if( temp[0] == "index" ) continue;

            boost::iter_split( temp, temp[0], boost::algorithm::first_finder( "-" ));
            samples.emplace( temp[0] );
        }

        for( auto& smp : samples )
        {
            if( boost::filesystem::exists( p + smp + ".tsv" )) boost::filesystem::remove( p + smp + ".tsv" );
            boost::filesystem::create_symlink(( "../LenDist/" + smp + "-isomiRs.tsv" ).c_str(), ( p + smp + ".tsv" ).c_str() );
        }
    }

    static void output_bubplot_visualization(
            const std::string& output_name,
            const std::string& node_path,
            const std::string& heatbub_js,
            const std::size_t& min_len,
            const std::size_t& max_len
            )
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
        output << "        $Chart_Types = $_POST['Chart_Types'];" << "\n";
        output << "        $GMPMT_Types = $_POST['GMPMT_Types'];" << "\n";
        output << "        $Annotation_Select = $_POST['Annotation_Select'];" << "\n";
        output << "        $Annotation_Arms = $_POST['Annotation_Arms'];" << "\n";
        output << "        $Sample_Files = $_POST['Sample_Files'];" << "\n";
        output << "        $Sequence_5p = $_POST['Sequence_5p'];" << "\n";
        output << "        $Sequence_3p = $_POST['Sequence_3p'];" << "\n";
        output << "        $Min_Length = $_POST['Min_Length'];" << "\n";
        output << "        $Max_Length = $_POST['Max_Length'];" << "\n";
        output << "        $isLog2 = $_POST['isLog2'];" << "\n";
        output << "        $Filter_Pref = $_POST['Filter_Pref'];" << "\n";
        output << "        $Filter_Type = $_POST['Filter_Type'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://code.jquery.com/jquery-3.3.1.min.js ></script>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== Filter ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Pref_List = array( 'UniLike', 'DisLike' );" << "\n";
        output << "        $Type_List = array( 'GMPM', 'GM', 'PM', 'Atail', 'Ctail', 'Gtail', 'Utail' );" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Filter_Pref onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option value=\"\"'; if($Filter_Pref=='') echo 'selected'; echo '>Preference Filter</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Pref_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Pref_List[$i].' ';" << "\n";
        output << "            if( $Filter_Pref == $Pref_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Pref_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo '</select>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=Filter_Type onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option value=\"\"'; if($Filter_Type=='') echo 'selected'; echo '>Preference Type</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Type_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Type_List[$i].' ';" << "\n";
        output << "            if( $Filter_Type == $Type_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Type_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Read Filter ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Filter_Anno = Array();" << "\n";
        output << "" << "\n";
        output << "        If( $Filter_Type != '' && $Filter_Type != 'Preference Type' &&" << "\n";
        output << "            $Filter_Pref != '' && $Filter_Pref != 'Preference Filter' )" << "\n";
        output << "        {" << "\n";
        output << "            $Samples = [ 'AGO1', 'AGO2', 'AGO3' ];" << "\n";
        output << "            For( $i = 0; $i < Count( $Samples ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $Filter_TSV = '../Preference/'.$Filter_Type.'/'.$Samples[ $i ].'_'.$Filter_Pref.'.text';" << "\n";
        output << "" << "\n";
        output << "                $inFile = new SplFileObject( $Filter_TSV );" << "\n";
        output << "                $Filter_Anno_AGO[ $i ] = Array();" << "\n";
        output << "" << "\n";
        output << "                while( !$inFile->eof() )" << "\n";
        output << "                {" << "\n";
        output << "                    $inFile_Lines = Rtrim( $inFile->fgets() );" << "\n";
        output << "                    if( $inFile_Lines == '' ) continue;" << "\n";
        output << "" << "\n";
        output << "                    $inFile_Lines = Explode( '*', $inFile_Lines )[0];" << "\n";
        output << "                    $inFile_Lines = Explode( '!', $inFile_Lines )[0];" << "\n";
        output << "                    $inFile_Lines = Explode( '|', $inFile_Lines )[0];" << "\n";
        output << "                    $inFile_Lines = Substr( $inFile_Lines, 0, StrLen( $inFile_Lines ) -3 );" << "\n";
        output << "" << "\n";
        output << "                    Array_Push( $Filter_Anno, $inFile_Lines );" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
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
        output << "            $anno = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "            $Anno_Array[ $anno[0] ] += $anno[1];" << "\n";
        output << "" << "\n";
        output << "            if( $Annotation_Select == $anno[0] )" << "\n";
        output << "            {" << "\n";
        output << "                for( $k = 2; $k < Count( $anno ); ++$k )" << "\n";
        output << "                    Array_Push( $Arm_Array, $anno[$k] );" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Total_PPM = 0.0;" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "            if( $Total_PPM < $ppm ) $Total_PPM = $ppm;" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Annotation_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $ppm )" << "\n";
        output << "        {" << "\n";
        output << "            If( $Filter_Type != '' && $Filter_Type != 'Preference Type'   &&" << "\n";
        output << "                $Filter_Pref != '' && $Filter_Pref != 'Preference Filter' &&" << "\n";
        output << "                !In_Array( $anno, $Filter_Anno )) continue;" << "\n";
        output << "" << "\n";
        output << "            echo '<option value='.$anno.' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Annotation_Select == $anno ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$anno.' ('.$ppm.'ppm)</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Chart_Types' value='$Chart_Types' />" << "\n";
        output << "            <input type='hidden' name='GMPMT_Types' value='$GMPMT_Types' />" << "\n";
        output << "            <input type='hidden' name='Sample_Files' value='$Sample_Files' />" << "\n";
        output << "            <input type='hidden' name='Filter_Pref' value='$Filter_Pref' />" << "\n";
        output << "            <input type='hidden' name='Filter_Type' value='$Filter_Type' />" << "\n";
        output << "            <input type='hidden' name='Min_Length' value='$Min_Length' />" << "\n";
        output << "            <input type='hidden' name='Max_Length' value='$Max_Length' />" << "\n";
        output << "            <input type='hidden' name='isLog2' value='$isLog2' />" << "\n";
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
        output << "#<!--================== ChartTypes Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Chart_Types onchange=this.form.submit();>';" << "\n";
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
        output << "        echo '</select>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPMT_Types Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        if( $Chart_Types == 'heatmap' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<select name=GMPMT_Types onchange=this.form.submit();>';" << "\n";
        output << "            echo \"<option value='' \"; if($GMPMT_Types=='') echo 'selected'; echo '>Types</option>';" << "\n";
        output << "" << "\n";
        output << "            $GMPMT_List = array('GMPM', 'GM', 'PM', 'A_Tail', 'C_Tail', 'G_Tail', 'U_Tail', 'Other_Tail');" << "\n";
        output << "" << "\n";
        output << "            For( $i = 0; $i < Count( $GMPMT_List ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value='.$GMPMT_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $GMPMT_Types == $GMPMT_List[$i] )" << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "                echo '>' . $GMPMT_List[$i] . '</option>';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '</select>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MinLength input =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text onchange=this.form.submit(); name=Min_Length size=3 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Min_Length == '' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '" << min_len << "';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $Min_Length;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== MaxLength input =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text onchange=this.form.submit(); name=Max_Length size=3 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Max_Length == '' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '" << max_len << "';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $Max_Length;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--================== AnnotatArm Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $Arms_List = array('5p', '3p');" << "\n";
        output << "        $Arms_Temp = array();" << "\n";
        output << "        $isSequBox = false;" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Arm_Array ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=checkbox onchange=this.form.submit(); name=Annotation_Arms[] value='.$Arms_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $Arm_Array[$i] != '' )" << "\n";
        output << "            {" << "\n";
        output << "                For( $j = 0; $j < Count( $Annotation_Arms ); ++$j )" << "\n";
        output << "                    if( $Annotation_Arms[$j] == $Arms_List[$i] )" << "\n";
        output << "                    {" << "\n";
        output << "                        Array_Push( $Arms_Temp, $Arms_List[$i] );" << "\n";
        output << "                        $isSequBox = true;" << "\n";
        output << "                        echo 'checked ';" << "\n";
        output << "                    }" << "\n";
        output << "            }" << "\n";
        output << "            else echo 'disabled=\"disabled\" ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$Arms_List[$i].'</input>';" << "\n";
        output << "" << "\n";
        output << "            if( $isSequBox )" << "\n";
        output << "            {" << "\n";
        output << "                $isSequBox = false;" << "\n";
        output << "                echo '<input type=text onchange=this.form.submit(); name=Sequence_'.$Arms_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "                if( $Arms_List[$i] == '5p' )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Sequence_5p == '' || $Sequence_5p == 'Sequence_5p' ) $Sequence_5p = $Arm_Array[$i];" << "\n";
        output << "                    echo 'size='.Strlen( $Sequence_5p )*1.4.' value='.$Sequence_5p.' />';" << "\n";
        output << "                }" << "\n";
        output << "                else" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Sequence_3p == '' || $Sequence_3p == 'Sequence_3p' ) $Sequence_3p = $Arm_Array[$i];" << "\n";
        output << "                    echo 'size='.Strlen( $Sequence_3p )*1.4.' value='.$Sequence_3p.' />';" << "\n";
        output << "                }" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $Annotation_Arms = $Arms_Temp;" << "\n";
        output << "" << "\n";
        output << "#<!--================== SampleFile Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv | grep -v AnnoSeq' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<input type=checkbox onchange=this.form.submit(); name=Sample_Files[] value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            For( $j = 0; $j < Count( $Sample_Files ); ++$j )" << "\n";
        output << "                if( $Sample_Files[$j] == $TSV_List[$i] )" << "\n";
        output << "                    echo 'checked ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.Explode( \".\", $TSV_List[$i] )[0].'</input>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--======================== isLog2 =========================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=checkbox onchange=this.form.submit(); name=isLog2 value=\"--log2\" ';" << "\n";
        output << "" << "\n";
        output << "        if( $isLog2 == \"--log2\" ) echo 'checked ';" << "\n";
        output << "        echo '>'.Log2.'</input>';" << "\n";
        output << "" << "\n";
        output << "        echo \"<input type='hidden' name='Annotation_Select' value='$Annotation_Select' /></form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--======================== Debug ==========================-->" << "\n";
        output << "" << "\n";
        output << "        $Script = '" << node_path << " " << heatbub_js << "';" << "\n";
        output << "        if( $Chart_Types != '' ) $Script = $Script.' --mode '.$Chart_Types;" << "\n";
        output << "        if( $GMPMT_Types != '' ) $Script = $Script.' --type '.$GMPMT_Types;" << "\n";
        output << "        if( $Annotation_Select != '' ) $Script = $Script.' --input '.$Annotation_Select;" << "\n";
        output << "        if( $Min_Length != '' && $Min_Length != 'minLen' ) $Script = $Script.' --minlen '.$Min_Length;" << "\n";
        output << "        if( $Max_Length != '' && $Max_Length != 'maxLen' ) $Script = $Script.' --maxlen '.$Max_Length;" << "\n";
        output << "        if( !Empty( $Annotation_Arms ))" << "\n"; 
        output << "        {" << "\n";
        output << "            $Script = $Script.' --arm '.Implode( '', $Annotation_Arms );" << "\n";
        output << "            For( $i = 0; $i < Count( $Annotation_Arms ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Annotation_Arms[ $i ] == '5p' ) $Script = $Script.' --arm5seq '.$Sequence_5p;" << "\n";
        output << "                if( $Annotation_Arms[ $i ] == '3p' ) $Script = $Script.' --arm3seq '.$Sequence_3p;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "        if( !Empty( $Sample_Files )) $Script = $Script.' --files '.Implode( ' ', $Sample_Files );" << "\n";
        output << "        if( $isLog2 != '' ) $Script = $Script.' '.$isLog2;" << "\n";
        output << "        echo '<br/>';" << "\n";
        output << "" << "\n";
        output << "        if(( $Chart_Types == 'bubble'  ||" << "\n";
        output << "           ( $Chart_Types == 'heatmap' && $GMPMT_Types != '' )) && $Annotation_Select != '' &&" << "\n";
        output << "             !Empty( $Annotation_Arms ) && !Empty( $Sample_Files ))" << "\n";
        output << "        {" << "\n";
        output << "            $Handle = Popen( \"$Script 2>&1\", 'r' );" << "\n";
        output << "            while( $Read = Fread( $Handle, 20096 )) $Response[] = Trim( $Read );" << "\n";
        output << "            Pclose( $Handle ); Flush();" << "\n";
        output << "            for( $i = 0; $i < Count( $Response ); ++$i ) Echo $Response[ $i ];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
