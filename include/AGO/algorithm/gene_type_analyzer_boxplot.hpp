#pragma once
#include <cstdlib>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBoxPlot
{

  public:

    GeneTypeAnalyzerBoxPlot()
    {}

    static void make_hete_table(
            const std::string& output_path,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::string& biotype,
            const std::string& token = ""
            )
    {
        std::vector< std::pair< 
            std::map< std::string, std::map< std::size_t, double >>, // 5p
            std::map< std::string, std::map< std::size_t, double >>  // 3p
        >> hete_tables( bed_samples.size() );

        std::size_t arm;
        std::string gene_name;

        bool is_arms = token == "3p" || token == "5p" ? true : false;
        bool is_lens = token != "" && !is_arms ? true : false;

        std::size_t* end_5p;
        std::size_t* end_3p;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                if( raw_bed.annotation_info_.empty() || raw_bed.annotation_info_[0].empty() ) continue;
                if( biotype == "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != "miRNA" && raw_bed.annotation_info_[0][0] != "mirtron" ) continue;
                if( biotype != "miRNA_mirtron" && raw_bed.annotation_info_[0][0] != biotype ) continue;

                for( int i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    gene_name = raw_bed.annotation_info_[0][ i+1 ];
                    arm = std::stoi( GeneTypeAnalyzerCounting::get_arm( gene_name ).substr( 0, 1 ));

                    if( is_arms && token != ( std::to_string( arm ) + "p"     )) continue;
                    if( is_lens && token !=   std::to_string( raw_bed.length_ )) continue;

                    switch( raw_bed.strand_ )
                    {
                        case '+' :
                            end_5p = &raw_bed.start_;
                            end_3p = &raw_bed.end_;
                            break;

                        case '-' :
                            end_5p = &raw_bed.end_;
                            end_3p = &raw_bed.start_;
                            break;
                    }

                    if( hete_tables[ smp ].first[ gene_name ].find( *end_5p ) == hete_tables[ smp ].first[ gene_name ].end() )
                        hete_tables[ smp ].first[ gene_name ][ *end_5p ] = 0;

                    if( hete_tables[ smp ].second[ gene_name ].find( *end_3p ) == hete_tables[ smp ].second[ gene_name ].end() )
                        hete_tables[ smp ].second[ gene_name ][ *end_3p ] = 0;

                    hete_tables[ smp ].first[ gene_name ][ *end_5p ] += raw_bed.ppm_;
                    hete_tables[ smp ].second[ gene_name ][ *end_3p ] += raw_bed.ppm_;
                }
            }

            format_hete_table( hete_tables[ smp ].first );
            format_hete_table( hete_tables[ smp ].second );
        }

        output_heterorgeneity( output_path, bed_samples, ano_len_idx, hete_tables, token );
    }

    static void format_hete_table( std::map< std::string, std::map< std::size_t, double >>& hete_tables )
    {
        double ppm_sum;
        std::size_t max_end;

        for( auto& anno : hete_tables )
        {
            max_end = get_max_end( anno.second );
            ppm_sum = get_ppm_sum( anno.second );

            for( auto& end : anno.second )
                end.second = (double)(std::abs( (int)(max_end) - (int)(end.first) )) * end.second;

            anno.second[0] = get_ppm_sum( anno.second ) / ppm_sum;
        }
    }

    static std::size_t get_max_end( std::map< std::size_t, double >& ends_table )
    {
        std::size_t max_end;
        double max_ppm = 0.0;

        for( auto& end : ends_table )
        {
            if( end.second > max_ppm )
            {
                max_end = end.first;
                max_ppm = end.second;
            }
        }

        return max_end;
    }

    static double get_ppm_sum( std::map< std::size_t, double >& ends_table )
    {
        double ppm_sum = 0.0;
        for( auto& end : ends_table ) ppm_sum += end.second;
        return ppm_sum;
    }

    static void output_heterorgeneity(
            const std::string& output_path,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            auto& hete_tables,
            const std::string& token
            )
    {
        std::vector< std::string > split;
        std::set< std::string > anno_idx;
        std::ofstream output;

        for( auto& anno : ano_len_idx.first )
        {
            boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));
            if( split.size() > 2 ) for( std::size_t i = 1; i < split.size() -2; ++i )
                    split[0] = split[0] + "_" + split[i];

            anno_idx.emplace( split[0] );
        }

        output.open( output_path + "Heterorgeneity_5p.tsv" );
        output << "Heterorgeneity_5p" << ( token == "" ? "" : ( "_" + token ));

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : anno_idx )
        {
            output << "\n" << anno;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( hete_tables[ smp ].first[ anno ].find( 0 ) != hete_tables[ smp ].first[ anno ].end() )
                {
                    output << "\t" << hete_tables[ smp ].first[ anno ][0];
                }
                else output << "\t0";
            }
        }

        output.close();

        output.open( output_path + "Heterorgeneity_3p.tsv" );
        output << "Heterorgeneity_3p" << ( token == "" ? "" : ( "_" + token ));

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : anno_idx )
        {
            output << "\n" << anno;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( hete_tables[ smp ].second[ anno ].find( 0 ) != hete_tables[ smp ].second[ anno ].end() )
                {
                    output << "\t" << hete_tables[ smp ].second[ anno ][0];
                }
                else output << "\t0";
            }
        }

        output.close();
    }

    static void make_rnafold_table(
            const std::string& output_path,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::string& rnafold_path,
            const std::string& biotype,
            const std::string& token = ""
            )
    {
        if(( biotype != "miRNA_mirtron" && biotype != "miRNA" && biotype != "mirtron" ) || rnafold_path == "" )
            return;

        std::set< std::string > annos;
        std::map< std::string, std::vector< std::string >> rnafolds;
        std::map< std::string, std::vector< std::map< std::string, double >>> entropies;
        std::map< std::string, std::vector< std::map< std::string, double >>> counts;

        GeneTypeAnalyzerSqalign::get_rnafold( rnafold_path, rnafolds );

        entropies[ "5p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "mid" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "3p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "3pTailOnly" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );

        counts[ "5p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "mid" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "3p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "3pTailOnly" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );

        bool is_arms = token == "3p" || token == "5p" ? true : false;
        bool is_lens = token != "" && !is_arms ? true : false;

        std::string arm;
        std::string mir;

        std::string gene_name;
        std::string gene_seed;

        std::size_t strpos;
        std::size_t midpos;
        std::size_t endpos;

        std::ofstream output;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                for( int i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    if( raw_bed.annotation_info_[0][i] != "miRNA" &&
                        raw_bed.annotation_info_[0][i] != "mirtron" &&
                        raw_bed.annotation_info_[0][i] != "miRNA_mirtron" )
                        continue;

                    gene_name = raw_bed.annotation_info_[0][ i+1 ];
                    arm = std::stoi( GeneTypeAnalyzerCounting::get_arm( gene_name ).substr( 0, 1 ));

                    gene_seed = raw_bed.getReadSeq( genome_table ).substr( 1, 7 )
                            + ( raw_bed.seed_md_tag != "" ? ( "|" + raw_bed.seed_md_tag ) : "" );

                    mir = gene_name.substr( 0, raw_bed.annotation_info_[0][1].length() -3 );
                    gene_name = gene_name + "_" + gene_seed;

                    if( is_arms && token != ( std::to_string( arm ) + "p"     )) continue;
                    if( is_lens && token !=   std::to_string( raw_bed.length_ )) continue;

                    if( raw_bed.start_ < std::stoi( rnafolds[ mir ][1] )) continue;
                    if( raw_bed.end_   > std::stoi( rnafolds[ mir ][2] )) continue;

                    if( entropies[ "5p"  ][ smp ].find( gene_name ) == entropies[ "5p"  ][ smp ].end() ) entropies[ "5p"  ][ smp ][ gene_name ] = 0.0;
                    if( entropies[ "mid" ][ smp ].find( gene_name ) == entropies[ "mid" ][ smp ].end() ) entropies[ "mid" ][ smp ][ gene_name ] = 0.0;
                    if( entropies[ "3p"  ][ smp ].find( gene_name ) == entropies[ "3p"  ][ smp ].end() ) entropies[ "3p"  ][ smp ][ gene_name ] = 0.0;

                    if( counts[ "5p"  ][ smp ].find( gene_name ) == counts[ "5p"  ][ smp ].end() ) counts[ "5p"  ][ smp ][ gene_name ] = 0.0;
                    if( counts[ "mid" ][ smp ].find( gene_name ) == counts[ "mid" ][ smp ].end() ) counts[ "mid" ][ smp ][ gene_name ] = 0.0;
                    if( counts[ "3p"  ][ smp ].find( gene_name ) == counts[ "3p"  ][ smp ].end() ) counts[ "3p"  ][ smp ][ gene_name ] = 0.0;

                    switch( raw_bed.strand_ )
                    {
                        case '+' : strpos = raw_bed.start_ - ( std::stoi( rnafolds[ mir ][1] ) + 1 ); break;
                        case '-' : strpos = ( std::stoi( rnafolds[ mir ][2] ) + 1 ) - raw_bed.end_  ; break;
                    }

                    annos.emplace( gene_name );
                    endpos = strpos + 1 + (int)raw_bed.length_ - (int)raw_bed.tail_length_;
                    std::vector< std::string > entropy( rnafolds[ mir ].begin() + 6, rnafolds[ mir ].begin() + rnafolds[ mir ].size() );

                    entropies[ "5p"  ][ smp ][ gene_name ] += get_entropy( strpos,     8, entropy );
                    entropies[ "mid" ][ smp ][ gene_name ] += get_entropy( strpos + 8, 4, entropy );
                    entropies[ "3p"  ][ smp ][ gene_name ] += get_entropy( endpos - 8, 8, entropy );

                    counts[ "5p"  ][ smp ][ gene_name ] += 1;
                    counts[ "mid" ][ smp ][ gene_name ] += 1;
                    counts[ "3p"  ][ smp ][ gene_name ] += 1;

                    if( raw_bed.getTail() != "" )
                    {
                        if( entropies[ "3pTailOnly" ][ smp ].find( gene_name ) == entropies[ "3pTailOnly" ][ smp ].end() ) entropies[ "3pTailOnly" ][ smp ][ gene_name ] = 0.0;
                        if( counts[ "3pTailOnly" ][ smp ].find( gene_name ) == counts[ "3pTailOnly" ][ smp ].end() ) counts[ "3pTailOnly" ][ smp ][ gene_name ] = 0.0;

                        entropies[ "3pTailOnly"  ][ smp ][ gene_name ] += get_entropy( endpos - 8, 8, entropy );
                        counts[ "3pTailOnly"  ][ smp ][ gene_name ] += 1;
                    }
                }
            }
        }

        for( auto& type : entropies )
        {
            output.open( output_path + "Entropy_" + type.first + ".tsv" );
            output << "Entropy_" << type.first;

            for( std::size_t smp = 0; smp < type.second.size(); ++smp )
                output << "\t" << bed_samples[ smp ].first;

            for( auto& anno : annos )
            {
                output << "\n" << anno;

                for( std::size_t smp = 0; smp < type.second.size(); ++smp )
                {
                    if( type.second[ smp ].find( anno ) != type.second[ smp ].end() )
                         output << "\t" << type.second[ smp ][ anno ] / counts[ type.first ][ smp ][ anno ];
                    else output << "\t0";
                }
            }

            output.close();
        }
    }

    static double get_entropy(
            std::size_t strpos,
            const std::size_t& length,
            const std::vector< std::string >& entropies
            )
    {
        double entropy = 0.0;

        if( strpos + length > entropies.size() )
            strpos = strpos - ( strpos + length - entropies.size() );

        for( std::size_t i = 0; i < length; i++ ) entropy += std::stod( entropies[ strpos + i ]);
        return entropy / (double)length;
    }

    static void debug( std::map< std::string, std::map< std::size_t, double >>& hete_table )
    {
        std::ofstream output( "./debug.tsv" );

        for( auto& gene : hete_table )
            for( auto& end : gene.second )
                output << gene.first << "\t" << end.first << "\t" << end.second << "\n";
        output.close();
    }

    static void output_boxplot_visualization( const std::string& output_name )
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
        output << "        $Type_List = array( 'GMPM', 'GM', 'PM', 'Atail', 'Ctail', 'Gtail', 'Ttail' );" << "\n";
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
