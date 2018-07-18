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
#include <AGO/algorithm/gene_type_analyzer_differential.hpp>
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
    std::string seedmap;
    std::string seedpie;
    std::string difference;
    std::string differential;

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
        , seedmap( "SeedMap/" )
        , seedpie( "SeedPie/" )
        , difference( "Difference/" )
        , differential( "Differential/" )
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
        , seedmap( "SeedMap/" )
        , seedpie( "SeedPie/" )
        , difference( "Difference/" )
        , differential( "Differential/" )
    {
        is_time_log = true;

        boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + taildot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + sa_plot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + volcano ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + barplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + ranking ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + difference ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + differential ));

        if( isSeed )
        {
            boost::filesystem::create_directory( boost::filesystem::path( output_path + seedmap ));
            boost::filesystem::create_directory( boost::filesystem::path( output_path + seedpie ));
        }
        else boost::filesystem::create_directory( boost::filesystem::path( output_path + sqalign ));

        std::chrono::time_point< std::chrono::system_clock > make_table_start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

        GeneTypeAnalyzerTaildot::make_taildot_table( biotype, ano_len_idx, bed_samples, genome_table, isSeed );
        GeneTypeAnalyzerSA_Plot::make_sa_plot_table( biotype, ano_len_idx, bed_samples, genome_table, isSeed );

        std::chrono::time_point< std::chrono::system_clock > make_table_end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
        if( is_time_log ) std::cerr << "make_table: " << std::chrono::duration< double >( make_table_end_time - make_table_start_time ).count() << "\n";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            const auto& sample_name = bed_samples[ smp ].first;

            if( !isSeed )
            {
                parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerDotplot::output_dotplot_isomirs( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "output_dotplot_isomirs " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    GeneTypeAnalyzerTaildot::output_taildot_isomirs( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "output_taildot_isomirs " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });
            }

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_dotplot " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerTaildot::output_taildot( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_taildot " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerSA_Plot::output_sa_plot( output_path + sa_plot, ano_len_idx, sample_name, smp );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_sa_plot " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerLendist::output_lendist( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_lendist " << smp << ": " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GMPM"    );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_barplot GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_barplot GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "PM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_barplot PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "Tailing" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_barplot Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_ranking GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_ranking GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_ranking PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_ranking Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });


        parallel_pool.job_post([ &biotype, &isSeed, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDotplot::output_dotplot_visualization( output_path + dotplot, biotype, isSeed );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_dotplot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &biotype, &isSeed, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerTaildot::output_taildot_visualization( output_path + taildot, biotype, isSeed );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_taildot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &biotype, &isSeed, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerSA_Plot::output_sa_plot_visualization( output_path + sa_plot );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_sa_plot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &isSeed, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist, isSeed );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_lendist_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerBarplot::output_barplot_visualization( output_path + barplot );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_barplot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + ranking );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_valplot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        if( !isSeed )
        {
            parallel_pool.job_post([ &bed_samples, &biotype, &genome_table, &extend_refseq, &max_anno_merge_size, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerSqalign::output_sqalign( output_path + sqalign, bed_samples, biotype, genome_table, extend_refseq, max_anno_merge_size );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_sqalign: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &biotype, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerMDTCpos::make_mdtcpos( bed_samples, biotype );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "make_mdtcpos: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &extend_refseq, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerSqalign::output_sqalign_visualization( output_path + sqalign, extend_refseq );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_sqalign_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerMDTCpos::output_mdtcpos_visualization( output_path );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_mdtcpos_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }
        else
        {
            parallel_pool.job_post([ &bed_samples, &seed_match_table, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerSeedMap::output_seed_mapping_table( output_path + seedmap, bed_samples, seed_match_table );
                GeneTypeAnalyzerSeedPie::output_seedpie_visualization( output_path + seedpie );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_seedpie_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerSeedMap::output_seedmap_visualization( output_path + seedmap );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_seedmap_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_loading_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_loading_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_loading_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_loading_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_length_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_length_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_length_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_length_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

            GeneTypeAnalyzerDifferential::output_loading_differential( output_path + differential, bed_samples, ano_len_idx, anno_table_tail );
            GeneTypeAnalyzerVolcano::output_volcano_visualization( output_path + volcano, biotype );

            std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
            if( is_time_log ) std::cerr << "output_volcano_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
        });



        if( biotype.substr( 0, 5 ) == "miRNA" || biotype == "mirtron" )
        {
            if( !isSeed )
            {
                boost::filesystem::create_directory( boost::filesystem::path( output_path + bubplot ));

                parallel_pool.job_post([ &bed_samples, &node_path, &heatbub_js, &min_len, &max_len, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    algorithm::GeneTypeAnalyzerBubplot::output_bubplot_visualization( output_path + bubplot, node_path, heatbub_js, min_len, max_len );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "output_bubplot_visualization: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });

                parallel_pool.job_post([ &bed_samples, &biotype, &thread_number, &extend_merge, &genome_table, this ] ()
                {
                    std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                    double ppm_filter = 1;
                    algorithm::GeneTypeAnalyzerBubplot::output_bubplot( output_path + bubplot, bed_samples, biotype, thread_number, extend_merge, ppm_filter, genome_table );

                    std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                    if( is_time_log ) std::cerr << "output_bubplot: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
                });
            }



            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_arms_difference GMPM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_arms_difference GM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_arms_difference PM: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                std::chrono::time_point< std::chrono::system_clock > start_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );

                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );

                std::chrono::time_point< std::chrono::system_clock > end_time = std::chrono::time_point< std::chrono::system_clock >( std::chrono::system_clock::now() );
                if( is_time_log ) std::cerr << "output_arms_difference Tailing: " << std::chrono::duration< double >( end_time - start_time ).count() << "\n";
            });
        }

        parallel_pool.flush_pool();
    }
};

} // end of namespace algorithm
} // end of namespace ago
