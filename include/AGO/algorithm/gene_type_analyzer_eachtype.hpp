#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_dotplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_taildot.hpp>
#include <AGO/algorithm/gene_type_analyzer_lendist.hpp>
#include <AGO/algorithm/gene_type_analyzer_barplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_bubplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_sqalign.hpp>
#include <AGO/algorithm/gene_type_analyzer_valplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_ranking.hpp>
#include <AGO/algorithm/gene_type_analyzer_mdtcpos.hpp>
#include <AGO/algorithm/gene_type_analyzer_difference.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerEachtype
{
    std::string output_path;
    ParaThreadPool parallel_pool;

    std::string dotplot;
    std::string taildot;
    std::string lendist;
    std::string barplot;
    std::string bubplot;
    std::string sqalign;
    std::string valplot;
    std::string ranking;
    std::string difference;

  public:

    GeneTypeAnalyzerEachtype()
        : output_path( "" )
        , parallel_pool( 0 )
        , dotplot( "DotPlot/" )
        , taildot( "TailDot/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , bubplot( "BubPlot/" )
        , sqalign( "SqAlign/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , difference( "Difference/" )
    {}

    GeneTypeAnalyzerEachtype(
            const std::string& biotype,
            std::string output_path_,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            std::size_t& thread_number,
            std::map< std::string, std::string >& genome_table,
            std::string& node_path,
            std::string& heatbub_js,
            std::size_t& min_len,
            std::size_t& max_len,
            std::size_t& extend_merge,
            std::size_t& extend_refseq
            )
        : output_path( output_path_
                + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" )
                + ( biotype == "miRNA_mirtron" ? "" : biotype )
                + "/" )
        , parallel_pool( thread_number )
        , dotplot( "DotPlot/" )
        , taildot( "TailDot/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , bubplot( "BubPlot/" )
        , sqalign( "SqAlign/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , difference( "Difference/" )
    {
        boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + taildot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + barplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + sqalign ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + ranking ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + difference ));

        GeneTypeAnalyzerTaildot::make_taildot_table( biotype, ano_len_idx, bed_samples, genome_table );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            const auto& sample_name = bed_samples[ smp ].first;

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerDotplot::output_dotplot_isomirs( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerTaildot::output_taildot_isomirs( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerTaildot::output_taildot( output_path + taildot, ano_len_idx, anno_mark[ smp ], sample_name, smp );
            });

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerLendist::output_lendist( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
            });
        }



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GMPM"    );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "GM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "PM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &biotype, this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot( output_path + barplot, bed_samples, ano_len_idx, anno_table_tail, biotype, "Tailing" );
        });



        parallel_pool.job_post([ &bed_samples, &biotype, &genome_table, &extend_refseq, this ] ()
        {
            GeneTypeAnalyzerSqalign::output_sqalign( output_path + sqalign, bed_samples, biotype, genome_table, extend_refseq );
        });

        parallel_pool.job_post([ &bed_samples, &biotype, this ] ()
        {
            GeneTypeAnalyzerMDTCpos::make_mdtcpos( bed_samples, biotype );
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
        {
            GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );
        });


        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerDotplot::output_dotplot_visualization( output_path + dotplot );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerTaildot::output_taildot_visualization( output_path + taildot );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot_visualization( output_path + barplot );
        });

        parallel_pool.job_post([ &extend_refseq, this ] ()
        {
            GeneTypeAnalyzerSqalign::output_sqalign_visualization( output_path + sqalign, extend_refseq );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + ranking );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerMDTCpos::output_mdtcpos_visualization( output_path );
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_loading_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );
        });



        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );
        });

        parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
        {
            GeneTypeAnalyzerDifference::output_length_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );
        });



        if( biotype.substr( 0, 5 ) == "miRNA" || biotype == "mirtron" )
        {
            boost::filesystem::create_directory( boost::filesystem::path( output_path + bubplot ));

            parallel_pool.job_post([ &bed_samples, &node_path, &heatbub_js, &min_len, &max_len, this ] ()
            {
                algorithm::GeneTypeAnalyzerBubplot::output_bubplot_visualization( output_path + bubplot, node_path, heatbub_js, min_len, max_len );
            });

            parallel_pool.job_post([ &bed_samples, &biotype, &thread_number, &extend_merge, &genome_table, this ] ()
            {
                double ppm_filter = 1;
                algorithm::GeneTypeAnalyzerBubplot::output_bubplot( output_path + bubplot, bed_samples, biotype, thread_number, extend_merge, ppm_filter, genome_table );
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GMPM" );
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "GM" );
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "PM" );
            });

            parallel_pool.job_post([ &bed_samples, &ano_len_idx, &anno_table_tail, this ] ()
            {
                GeneTypeAnalyzerDifference::output_arms_difference( output_path + difference, bed_samples, ano_len_idx, anno_table_tail, "Tailing" );
            });
        }

        parallel_pool.flush_pool();
    }
};

} // end of namespace algorithm
} // end of namespace ago
