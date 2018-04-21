#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_dotplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_lendist.hpp>
#include <AGO/algorithm/gene_type_analyzer_barplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_valplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_ranking.hpp>
#include <AGO/algorithm/gene_type_analyzer_difference.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerEachtype
{
    std::string output_path;
    ParaThreadPool parallel_pool;

    std::string dotplot;
    std::string lendist;
    std::string barplot;
    std::string valplot;
    std::string ranking;
    std::string difference;

  public:

    GeneTypeAnalyzerEachtype()
        : output_path( "" )
        , parallel_pool( 0 )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , difference( "Difference/" )
    {}

    GeneTypeAnalyzerEachtype(
            std::string& biotype,
            std::string output_path_,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark,
            std::size_t& thread_number
            )
        : output_path( output_path_ + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" ) + biotype + "/" )
        , parallel_pool( thread_number )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , barplot( "BarPlot/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , difference( "Difference/" )
    {
        boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + barplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + ranking ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + difference ));

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            const auto& sample_name = bed_samples[ smp ].first;

            parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
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
            GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerBarplot::output_barplot_visualization( output_path + barplot );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot );
        });

        parallel_pool.job_post([ this ] ()
        {
            GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + ranking );
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



        if( biotype == "miRNA" )
        {
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
