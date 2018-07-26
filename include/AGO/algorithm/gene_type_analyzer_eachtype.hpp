#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_dotplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_taildot.hpp>
#include <AGO/algorithm/gene_type_analyzer_sa_plot.hpp>
#include <AGO/algorithm/gene_type_analyzer_volcano.hpp>
#include <AGO/algorithm/gene_type_analyzer_lendist.hpp>
#include <AGO/algorithm/gene_type_analyzer_barplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_bubplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_sqalign.hpp>
#include <AGO/algorithm/gene_type_analyzer_valplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_ranking.hpp>
#include <AGO/algorithm/gene_type_analyzer_seedmap.hpp>
#include <AGO/algorithm/gene_type_analyzer_seedpie.hpp>
#include <AGO/algorithm/gene_type_analyzer_mdtcpos.hpp>
#include <AGO/algorithm/gene_type_analyzer_difference.hpp>
#include <AGO/algorithm/gene_type_analyzer_diffbar.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>
#include <chrono>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerEachtype
{
    bool is_time_log;
    std::string output_path;
    ParaThreadPool parallel_pool;

    std::string dotplot;
    std::string taildot;
    std::string sa_plot;
    std::string volcano;
    std::string lendist;
    std::string barplot;
    std::string bubplot;
    std::string sqalign;
    std::string valplot;
    std::string ranking;
    std::string mdtcpos;
    std::string seedmap;
    std::string seedpie;
    std::string diffbar;
    std::string difference;

  public:

    GeneTypeAnalyzerEachtype()
        : is_time_log( false )
        , output_path( "" )
        , parallel_pool( 0 )
        , dotplot( "DotPlot/" )
        , taildot( "TailDot/" )
        , sa_plot( "SA_Plot/" )
        , volcano( "Volcano/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , bubplot( "BubPlot/" )
        , sqalign( "SqAlign/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , mdtcpos( "MDTCpos/" )
        , seedmap( "SeedMap/" )
        , seedpie( "SeedPie/" )
        , diffbar( "DiffBar/" )
        , difference( "Difference/" )
    {}

    GeneTypeAnalyzerEachtype(
            const std::string& biotype,
            const std::string output_path_,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::map< std::string, double >>>& seed_match_table,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            std::map< std::string, std::string >& genome_table,
            const std::size_t& thread_number,
            const std::string& node_path,
            const std::string& heatbub_js,
            const std::size_t& min_len,
            const std::size_t& max_len,
            const std::size_t& extend_merge,
            const std::size_t& extend_refseq,
            const std::size_t& max_anno_merge_size,
            const bool& webpage_update_only,
            const bool& isSeed
            )
        : is_time_log( false )
        , output_path( output_path_
                + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" )
                + ( biotype == "miRNA_mirtron" ? "" : biotype )
                + "/" )
        , parallel_pool( thread_number )
        , dotplot( "DotPlot/" )
        , taildot( "TailDot/" )
        , sa_plot( "SA_Plot/" )
        , volcano( "Volcano/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , bubplot( "BubPlot/" )
        , sqalign( "SqAlign/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , mdtcpos( "MDTCpos/" )
        , seedmap( "SeedMap/" )
        , seedpie( "SeedPie/" )
        , diffbar( "DiffBar/" )
        , difference( "Difference/" )
    {
        if( !webpage_update_only )
        {
            // is_time_log = true;

            boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + taildot ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + sa_plot ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + volcano ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + ranking ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + diffbar ));
            // boost::filesystem::create_directory( boost::filesystem::path( output_path + barplot ));
            // boost::filesystem::create_directory( boost::filesystem::path( output_path + difference ));

            if( isSeed )
            {
                boost::filesystem::create_directory( boost::filesystem::path( output_path + seedmap ));
                boost::filesystem::create_directory( boost::filesystem::path( output_path + seedpie ));
            }
            else
            {
                boost::filesystem::create_directory( boost::filesystem::path( output_path + mdtcpos ));
                boost::filesystem::create_directory( boost::filesystem::path( output_path + sqalign ));
            }

            std::chrono::time_point< std::chrono::system_clock > make_table_start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerTaildot::make_taildot_table( biotype, ano_len_idx, bed_samples, genome_table, isSeed );
            GeneTypeAnalyzerSA_Plot::make_sa_plot_table( biotype, ano_len_idx, bed_samples, genome_table, isSeed );

            std::chrono::time_point< std::chrono::system_clock > make_table_end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "make_table: " << std::chrono::duration< double >( make_table_end_time - make_table_start_time ).count() << "\n";

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                const auto& sample_name = bed_samples[ smp ].first;

                parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, &isSeed, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    if( !isSeed )
                    {
                        GeneTypeAnalyzerDotplot::output_dotplot_isomirs( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
                        GeneTypeAnalyzerTaildot::output_taildot_isomirs( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );
                        GeneTypeAnalyzerLendist::output_lendist_isomirs( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
                    }

                    GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
                    GeneTypeAnalyzerTaildot::output_taildot( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );
                    GeneTypeAnalyzerLendist::output_lendist( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "LenDist & DotPlot & TailDot" << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerSA_Plot::output_sa_plot( output_path + sa_plot, ano_len_idx, sample_name, smp );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "SA_Plot " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });
            }

            // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
            // {
            //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            //     GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GMPM"    );
            //     GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GM"      );
            //     GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "PM"      );
            //     GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "Tailing" );

            //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            //     if( is_time_log ) std::cerr << "BarPlot: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            // });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "GMPM"    );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "GM"      );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "PM"      );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Tailing" );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Atail"   );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Ctail"   );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Gtail"   );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Ttail"   );
                GeneTypeAnalyzerDiffBar::output_diffbar( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Other"   );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "DiffBar: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
                GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
                GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
                GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "ValPlot: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
                GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
                GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
                GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "Ranking : " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDiffBar::output_loading_differential( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, false );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "Differential: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }

        if( !isSeed )
        {
            if( !webpage_update_only )
            {
                parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "GMPM"    );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "GM"      );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "PM"      );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Tailing" );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Atail"   );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Ctail"   );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Gtail"   );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Ttail"   );
                    GeneTypeAnalyzerDiffBar::output_diffbar_isomirs( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, biotype, "Other"   );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "DiffBar_isomirs: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerValplot::output_valplot_isomirs( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
                    GeneTypeAnalyzerValplot::output_valplot_isomirs( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
                    GeneTypeAnalyzerValplot::output_valplot_isomirs( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
                    GeneTypeAnalyzerValplot::output_valplot_isomirs( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "ValPlot_isomirs: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerRanking::output_ranking_isomirs( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
                    GeneTypeAnalyzerRanking::output_ranking_isomirs( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
                    GeneTypeAnalyzerRanking::output_ranking_isomirs( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
                    GeneTypeAnalyzerRanking::output_ranking_isomirs( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "Ranking_isomirs : " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerDiffBar::output_loading_differential( output_path + diffbar, bed_samples, ano_len_idx, anno_table_tail, true );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "Differential_isomirs: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });
            }

            parallel_pool.job_post([ &bed_samples, &biotype, &genome_table, &extend_refseq, &max_anno_merge_size, &webpage_update_only, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                if( !webpage_update_only )
                    GeneTypeAnalyzerSqalign::output_sqalign( output_path + sqalign, bed_samples, biotype, genome_table, extend_refseq, max_anno_merge_size );

                GeneTypeAnalyzerSqalign::output_sqalign_visualization( output_path + sqalign, extend_refseq );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "SqAlign: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &biotype, &webpage_update_only, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerMDTCpos::make_mdtcpos( bed_samples, biotype );
                GeneTypeAnalyzerMDTCpos::output_mdtcpos_visualization( output_path + mdtcpos );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "MDTCpos: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }
        else
        {
            parallel_pool.job_post([ &bed_samples, &seed_match_table, &webpage_update_only, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                if( !webpage_update_only )
                    GeneTypeAnalyzerSeedMap::output_seed_mapping_table( output_path + seedmap, bed_samples, seed_match_table );

                GeneTypeAnalyzerSeedPie::output_seedpie_visualization( output_path + seedpie );
                GeneTypeAnalyzerSeedMap::output_seedmap_visualization( output_path + seedmap );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "SeedPie & SeedMap: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }

        parallel_pool.job_post([ &biotype, &isSeed, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDotplot::output_dotplot_visualization( output_path + dotplot, biotype, isSeed );
            GeneTypeAnalyzerTaildot::output_taildot_visualization( output_path + taildot, biotype, isSeed );
            GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist, isSeed );
            GeneTypeAnalyzerSA_Plot::output_sa_plot_visualization( output_path + sa_plot );
            // GeneTypeAnalyzerBarplot::output_barplot_visualization( output_path + barplot );
            GeneTypeAnalyzerDiffBar::output_diffbar_visualization( output_path + diffbar, isSeed );
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot, isSeed );
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + ranking, isSeed );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "Visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_loading_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_loading_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_loading_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_loading_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });



        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_length_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_length_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_length_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        // {
        //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        //     GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

        //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        //     if( is_time_log ) std::cerr << "output_length_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        // });

        if( biotype.substr( 0, 5 ) == "miRNA" || biotype == "mirtron" )
        {
            if( !isSeed )
            {
                if( !webpage_update_only )
                    boost::filesystem::create_directory( boost::filesystem::path( output_path + bubplot ));

                parallel_pool.job_post([ &bed_samples, &node_path, &heatbub_js, &min_len, &max_len, &biotype, &thread_number, &extend_merge, &genome_table, &webpage_update_only, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    double ppm_filter = 1;

                    if( !webpage_update_only )
                        algorithm::GeneTypeAnalyzerBubplot::output_bubplot( output_path + bubplot, bed_samples, biotype, thread_number, extend_merge, ppm_filter, genome_table );

                    algorithm::GeneTypeAnalyzerBubplot::output_bubplot_visualization( output_path + bubplot, node_path, heatbub_js, min_len, max_len );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "Bubplot: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });
            }

            // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            // {
            //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            //     GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

            //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            //     if( is_time_log ) std::cerr << "output_arms_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            // });

            // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            // {
            //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            //     GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

            //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            //     if( is_time_log ) std::cerr << "output_arms_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            // });

            // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            // {
            //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            //     GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

            //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            //     if( is_time_log ) std::cerr << "output_arms_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            // });

            // parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            // {
            //     std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            //     GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

            //     std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            //     if( is_time_log ) std::cerr << "output_arms_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            // });
        }

        parallel_pool.flush_pool();

        GeneTypeAnalyzerVolcano::output_volcano_visualization( output_path + volcano, biotype, isSeed );
    }
};

} // end of namespace algorithm
} // end of namespace ago
