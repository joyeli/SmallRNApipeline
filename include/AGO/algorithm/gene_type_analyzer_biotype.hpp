#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_dotplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_lendist.hpp>
#include <AGO/algorithm/gene_type_analyzer_valplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBiotype
{
    std::string output_path;
    ParaThreadPool smp_parallel_pool;
    AnnoLengthIndexType ano_len_idx;

    std::vector< std::vector< CountingTableType >> anno_table_tail;
    std::vector< std::map< std::string, std::string >> anno_mark;

    std::string biotype;
    std::string dotplot;
    std::string lendist;
    std::string valplot;

  public:

    GeneTypeAnalyzerBiotype()
        : output_path( "" )
        , smp_parallel_pool( 0 )
        , ano_len_idx( std::make_pair( std::set< std::string >(), std::set< std::size_t >() ))
        , anno_table_tail( 0, std::vector< CountingTableType >( 6, CountingTableType() ))
        , anno_mark( 1, std::map< std::string, std::string >() )
        , biotype( "BioType/" )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , valplot( "ValPlot/" )
    {}

    GeneTypeAnalyzerBiotype(
            const std::string output_path_,
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            std::vector< std::string >& biotype_list,
            const std::size_t& min_len,
            const std::size_t& max_len,
            const double& sudo_count
            )
        : output_path( output_path_ + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" ))
        , smp_parallel_pool( bed_samples.size() )
        , ano_len_idx( GeneTypeAnalyzerCounting::get_ano_len_idx( genome_table, bed_samples ))
        , anno_table_tail( bed_samples.size(), std::vector< CountingTableType >( 6, CountingTableType() ))
        , anno_mark( 1, std::map< std::string, std::string >() )
        , biotype( "BioType/" )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , valplot( "ValPlot/" )
    {
        boost::filesystem::create_directory( boost::filesystem::path( output_path + biotype ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &bed_samples, &genome_table, this ] ()
            {
                GeneTypeAnalyzerCounting::make_anno_table( bed_samples[ smp ].second, anno_table_tail[ smp ], anno_mark[0], genome_table );
            });
        }

        smp_parallel_pool.flush_pool();

        for( auto& biotype : biotype_list )
        {
            if( biotype == "miRNA_mirtron" )
            {
                ano_len_idx.first.emplace( "miRNA" );
                ano_len_idx.first.emplace( "mirtron" );
            }
            else ano_len_idx.first.emplace( biotype );
        }

        GeneTypeAnalyzerCounting::table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            const auto& sample_name = bed_samples[ smp ].first;

            smp_parallel_pool.job_post([ smp, &sample_name, this ] ()
            {
                output_biotype( output_path + biotype + sample_name, ano_len_idx, anno_table_tail[ smp ], "GMPM" );
            });

            smp_parallel_pool.job_post([ smp, &sample_name, this ] ()
            {
                output_biotype( output_path + biotype + sample_name, ano_len_idx, anno_table_tail[ smp ], "GM" );
            });

            smp_parallel_pool.job_post([ smp, &sample_name, this ] ()
            {
                output_biotype( output_path + biotype + sample_name, ano_len_idx, anno_table_tail[ smp ], "PM" );
            });

            smp_parallel_pool.job_post([ smp, &sample_name, this ] ()
            {
                GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[0], sample_name );
            });

            smp_parallel_pool.job_post([ smp, &sample_name, this ] ()
            {
                GeneTypeAnalyzerLendist::output_lendist( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[0], sample_name );
            });
        }

        output_biotype_visualization( output_path + biotype );
        GeneTypeAnalyzerDotplot::output_dotplot_visualization( output_path + dotplot );
        GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist );

        output_biotype( output_path + biotype, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM" );
        output_biotype( output_path + biotype, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"   );
        output_biotype( output_path + biotype, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"   );

        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

        GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot );

        smp_parallel_pool.flush_pool();
    }

    void output_biotype(
            const std::string& output_name,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail, // 0-A / 1-C / 2-G / 3-T / 4-O / 5-GM
            const std::string& token
            )
    {
        std::ofstream output( output_name + "_" + token + ".tsv" );
        output << "Annotation";

        std::vector< double > gm( ano_len_idx.second.size(), 0.0 );
        std::vector< double > pm( ano_len_idx.second.size(), 0.0 );

        int ldx = -1;

        for( auto& len  : ano_len_idx.second ) output << "\t" << len;
        for( auto& anno : ano_len_idx.first )
        {
            output << "\n" << anno << "\t" << std::setprecision( 0 ) << std::fixed;

            gm = std::vector< double >( ano_len_idx.second.size(), 0.0 ); 
            pm = std::vector< double >( ano_len_idx.second.size(), 0.0 ); 

            if( token != "PM" )
            {
                if( anno_table_tail[5].find( anno ) != anno_table_tail[5].end() )
                {
                    ldx = -1;
                    for( auto& len : ano_len_idx.second )
                    {
                        ldx++;
                        if( anno_table_tail[5][ anno ].find( len ) == anno_table_tail[5][ anno ].end() ) continue;
                        gm[ ldx ] += anno_table_tail[5][ anno ][ len ];
                    }
                }
            }

            if( token != "GM" )
            {
                for( std::size_t i = 0; i < 5; i++ )
                {
                    if( anno_table_tail[i].find( anno ) != anno_table_tail[i].end() )
                    {
                        ldx = -1;
                        for( auto& len : ano_len_idx.second )
                        {
                            ldx++;
                            if( anno_table_tail[i][ anno ].find( len ) == anno_table_tail[i][ anno ].end() ) continue;
                            pm[ ldx ] += anno_table_tail[i][ anno ][ len ];
                        }
                    }
                }
            }

            for( std::size_t ldx = 0; ldx < ano_len_idx.second.size(); ++ldx )
                output << "\t" << ( token == "GMPM" ? gm[ ldx ] + pm[ ldx ] : ( token == "GM" ? gm[ ldx ] : pm[ ldx ] ));
        }

        output << "\n";
        output.close();
    }

    std::vector< double > get_sums( const std::vector< std::pair< std::string, std::vector< double >>>& value_vecs )
    {
        std::vector< double > sample_sums( value_vecs[0].second.size(), 0 );
        for( auto& anno : value_vecs )
        {
            for( std::size_t smp = 0; smp < anno.second.size(); ++smp )
                sample_sums[ smp ] += anno.second[ smp ];
        }
        return sample_sums;
    }

    void output_biotype(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            const std::string& token
            )
    {
        std::ofstream output( output_name + "All_" + token + ".tsv" );
        output << "Annotation";

        double gm = 0.0;
        double pm = 0.0;

        std::vector< double > sample_sums;
        std::vector< double > sample_values;
        std::vector< std::pair< std::string, std::vector< double >>> value_vecs;

        for( auto& smp  : bed_samples ) output << "\t" << smp.first;
        for( auto& anno : ano_len_idx.first )
        {
            sample_values = std::vector< double >( bed_samples.size() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                gm = 0.0;
                pm = 0.0;

                if( token != "PM" )
                {
                    if( anno_table_tail[ smp ][5].find( anno ) != anno_table_tail[ smp ][5].end() )
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table_tail[ smp ][5][ anno ].find( len ) != anno_table_tail[ smp ][5][ anno ].end() )
                                gm += anno_table_tail[ smp ][5][ anno ][ len ];
                        }
                }

                if( token != "GM" )
                {
                    for( std::size_t i = 0; i < 5; i++ )
                    {
                        if( anno_table_tail[ smp ][i].find( anno ) != anno_table_tail[ smp ][i].end() )
                            for( auto& len : ano_len_idx.second )
                            {
                                if( anno_table_tail[ smp ][i][ anno ].find( len ) != anno_table_tail[ smp ][i][ anno ].end() )
                                    pm += anno_table_tail[ smp ][i][ anno ][ len ];
                            }
                    }
                }

                sample_values[ smp ] = 
                    ( token == "GMPM" ? gm + pm : ( token == "GM" ? gm : ( token == "PM" ? pm : (( gm + pm ) < 1 ? 0 : ( pm * 100 / ( gm + pm ))))));
            }

            value_vecs.emplace_back( std::make_pair( anno +
                ( anno_mark[0].find( anno ) != anno_mark[0].end() ? anno_mark[0][ anno ] : "" ), sample_values ));
        }

        sample_sums = get_sums( value_vecs );

        for( auto& anno : value_vecs )
        {
            output << "\n" << anno.first << std::setprecision( 0 ) << std::fixed;
            for( std::size_t smp = 0; smp < anno.second.size(); ++smp )
                output << "\t" << anno.second[ smp ] * 100 / sample_sums[ smp ];
        }

        output << "\n";
        output.close();
    }

    void output_biotype_visualization( const std::string& output_name )
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
        output << "        $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "" << "\n";
        output << "        echo '<script src=https://d3js.org/d3.v3.js></script>';" << "\n";
        output << "        echo '<script src=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.min.js></script>';" << "\n";
        output << "        echo '<link href=https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.css rel=stylesheet type=text/css>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== GMPM ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=GMPM onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($GMPM=='') echo 'selected'; echo '>GM or PM</option>';" << "\n";
        output << "" << "\n";
        output << "        $GMPM_List = array('GMPM', 'GM', 'PM');" << "\n";
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
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "" << "\n";
        output << "#<!--================== Single TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep _'.$GMPM.'.tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='GMPM' value='$GMPM' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "        echo '<style>';" << "\n";
        output << "        echo 'text[ style=\"opacity: 0; text-anchor: middle;\" ]{';" << "\n";
        output << "        echo 'opacity: 1 !important';" << "\n";
        output << "        echo '}';" << "\n";
        output << "        echo '</style>';" << "\n";
        output << "" << "\n";
        output << "#<!--================== BioType ====================-->" << "\n";
        output << "" << "\n";
        output << "        $File_Tag = Explode( '_', $TSV_File );" << "\n";
        output << "" << "\n";
        output << "        echo \"<svg id='bar'></svg>" << "\n";
        output << "            <script>" << "\n";
        output << "                d3.tsv( '$TSV_File', function( tsv_data ) {" << "\n";
        output << "                    var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'Annotation'; });" << "\n";
        output << "                    var data = new Array();" << "\n";
        output << "                    tsv_data.forEach( function(d) {" << "\n";
        output << "                            d.len = length.map( function(key) {" << "\n";
        output << "                                return {" << "\n";
        output << "                                    x: key," << "\n";
        output << "                                    y: +d[key]" << "\n";
        output << "                                };" << "\n";
        output << "            " << "\n";
        output << "                            });" << "\n";
        output << "                            var sample = new Object();" << "\n";
        output << "                            sample.key = d['Annotation'];" << "\n";
        output << "                            sample.values = d.len;" << "\n";
        output << "                            data.push( sample );" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                nv.addGraph({" << "\n";
        output << "                        generate: function() {" << "\n";
        output << "                            var width = nv.utils.windowSize().width," << "\n";
        output << "                                height = nv.utils.windowSize().height;" << "\n";
        output << "                " << "\n";
        output << "                            var chart = nv.models.multiBarChart()" << "\n";
        output << "                                .margin({'left':70})" << "\n";
        output << "                                .width(width)" << "\n";
        output << "                                .height(height)" << "\n";
        output << "                                .stacked(true);" << "\n";
        output << "                " << "\n";
        output << "                            chart.dispatch.on('renderEnd', function(){" << "\n";
        output << "                                console.log('Render Complete');" << "\n";
        output << "                            });" << "\n";
        output << "            " << "\n";
        output << "                            var svg = d3.select('#bar').attr('height',height*2).datum(data);" << "\n";
        output << "                            console.log('calling chart');" << "\n";
        output << "                            svg.transition().duration(0).call(chart);" << "\n";
        output << "                " << "\n";
        output << "                            return chart;" << "\n";
        output << "                        }," << "\n";
        output << "                        callback: function(graph) {" << "\n";
        output << "                            nv.utils.windowResize(function() {" << "\n";
        output << "                                var width = nv.utils.windowSize().width;" << "\n";
        output << "                                var height = nv.utils.windowSize().height;" << "\n";
        output << "                                graph.width(width).height(height);" << "\n";
        output << "                " << "\n";
        output << "                                d3.select('#bar')" << "\n";
        output << "                                    .attr('width', width)" << "\n";
        output << "                                    .attr('height', height)" << "\n";
        output << "                                    .transition().duration(0)" << "\n";
        output << "                                    .call(graph);" << "\n";
        output << "                " << "\n";
        output << "                            });" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "                });" << "\n";
        output << "            </script>\";" << "\n";
        output << "        " << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";
        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
