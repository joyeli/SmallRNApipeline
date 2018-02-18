#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_dotplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_lendist.hpp>
#include <AGO/algorithm/gene_type_analyzer_valplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_ranking.hpp>
#include <AGO/algorithm/gene_type_analyzer_bubplot.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerEachtype
{
    std::string output_path;
    ParaThreadPool smp_parallel_pool;

    std::string dotplot;
    std::string lendist;
    std::string valplot;
    std::string ranking;
    std::string bubplot;

  public:

    GeneTypeAnalyzerEachtype()
        : output_path( "" )
        , smp_parallel_pool( 0 )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , bubplot( "Bubplot/" )
    {}

    GeneTypeAnalyzerEachtype(
            std::string output_path_,
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< std::map< std::string, std::string >>& anno_mark
            )
        : output_path( output_path_ + ( output_path_.at( output_path_.length() -1 ) != '/' ? "/" : "" ))
        , smp_parallel_pool( bed_samples.size() )
        , dotplot( "DotPlot/" )
        , lendist( "LenDist/" )
        , valplot( "ValPlot/" )
        , ranking( "Ranking/" )
        , bubplot( "Bubplot/" )
    {
        boost::filesystem::create_directory( boost::filesystem::path( output_path + dotplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + lendist ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + valplot ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + ranking ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + bubplot ));

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            const auto& sample_name = bed_samples[ smp ].first;

            smp_parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerDotplot::output_dotplot( output_path + dotplot, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
            });

            smp_parallel_pool.job_post([ smp, &sample_name, &ano_len_idx, &anno_table_tail, &anno_mark, this ] ()
            {
                GeneTypeAnalyzerLendist::output_lendist( output_path + lendist, ano_len_idx, anno_table_tail[ smp ], anno_mark[ smp ], sample_name );
            });
        }


        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
        GeneTypeAnalyzerValplot::output_valplot( output_path + valplot, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

        GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GMPM"    );
        GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "GM"      );
        GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "PM"      );
        GeneTypeAnalyzerRanking::output_ranking( output_path + ranking, bed_samples, ano_len_idx, anno_table_tail, anno_mark, "Tailing" );

        GeneTypeAnalyzerBubplot::output_bubplot( output_path + bubplot, bed_samples );

        GeneTypeAnalyzerDotplot::output_dotplot_visualization( output_path + dotplot );
        GeneTypeAnalyzerLendist::output_lendist_visualization( output_path + lendist );

        GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + valplot );
        GeneTypeAnalyzerValplot::output_valplot_visualization( output_path + ranking );

        smp_parallel_pool.flush_pool();
    }
};

} // end of namespace algorithm
} // end of namespace ago
