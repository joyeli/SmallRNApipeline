#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_filtering.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_biotype.hpp>
#include <AGO/algorithm/gene_type_analyzer_quantile.hpp>
#include <AGO/algorithm/gene_type_analyzer_eachtype.hpp>
#include <AGO/algorithm/gene_type_analyzer_annobed.hpp>
#include <AGO/algorithm/gene_type_analyzer_debug.hpp>
#include <AGO/algorithm/gene_type_analyzer_tmm.hpp>
#include <iomanip>

namespace ago {
namespace component {

class GeneTypeAnalyzer
    : public engine::NamedComponent
      , algorithm::GeneTypeAnalyzerFiltering
      , algorithm::GeneTypeAnalyzerCounting
      , algorithm::GeneTypeAnalyzerBiotype
      , algorithm::GeneTypeAnalyzerQuantile
      , algorithm::GeneTypeAnalyzerEachtype
      , algorithm::GeneTypeAnalyzerAnnobed
      , algorithm::GeneTypeAnalyzerTmm
{
    using Base = engine::NamedComponent;

    std::vector< std::string > biotype_list;
    std::vector< std::string > analysis_list;

    double sudo_count;

    bool output_annobed;
    bool is_skip_un_annotated;
    bool is_keep_other_biotype;
    bool webpage_update_only;

    std::size_t thread_number;
    std::size_t extend_refseq;
    std::size_t extend_merge;
    std::size_t max_anno_merge_size;
    std::size_t min_len;
    std::size_t max_len;

    std::string node_path;
    std::string heatbub_js;

    void emplace_list( const bpt::ptree& p, const std::string& tag, std::vector< std::string >& list )
    {
        std::size_t i = 0;
        std::size_t miRNA_idx = 0;
        std::size_t mirtron_idx = 0;
        std::size_t miRNA_mirtron_idx = 0;

        for( auto& biotype : p.get_child( tag ))
        {
            i++;
            if( biotype.second.data() == "miRNA" ) miRNA_idx = i; 
            if( biotype.second.data() == "mirtron" ) mirtron_idx = i; 
            if( biotype.second.data() == "miRNA_mirtron" ) miRNA_mirtron_idx = i; 
            list.emplace_back( biotype.second.data() );
        }

        if(( miRNA_mirtron_idx != 0 && ( miRNA_idx == 0 || mirtron_idx == 0 ))
        || ( miRNA_idx != 0 && mirtron_idx != 0 && miRNA_idx < mirtron_idx  ))
        {
            std::vector< std::string > biotype_list_temp;

            for( std::size_t i = 0; i < list.size(); ++i )
            {
                if( list[i] == "mirtron" ) continue;
                if( miRNA_mirtron_idx != 0 )
                {
                    if( list[i] == "miRNA" ) continue;
                    if( list[i] == "miRNA_mirtron" )
                    {
                        biotype_list_temp.emplace_back( "mirtron" );
                        biotype_list_temp.emplace_back( "miRNA" );
                    }
                }
                else if( list[i] == "miRNA" )
                {
                    biotype_list_temp.emplace_back( "mirtron" );
                    biotype_list_temp.emplace_back( "miRNA" );
                }

                biotype_list_temp.emplace_back( list[i] );
            }

            list = biotype_list_temp;
        }

        if( !is_skip_un_annotated )
            list.emplace_back( "un_annotated" );
    }

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        sudo_count = p.get_optional< double   >( "sudo_count"     ).value_or( 0.000001 );
        output_annobed = p.get_optional< bool >( "output_annobed" ).value_or( true     );
        is_skip_un_annotated  = p.get_optional< bool >( "is_skip_un_annotated"  ).value_or( false );
        is_keep_other_biotype = p.get_optional< bool >( "is_keep_other_biotype" ).value_or( false );
        webpage_update_only   = p.get_optional< bool >( "webpage_update_only"   ).value_or( false );
        thread_number  = p.get_optional< std::size_t >( "thread_number" ).value_or( 8  );
        extend_refseq  = p.get_optional< std::size_t >( "extend_refseq" ).value_or( 10 );
        extend_merge   = p.get_optional< std::size_t >( "extend_merge"  ).value_or( 2  );
        max_anno_merge_size = p.get_optional< std::size_t >( "max_anno_merge_size" ).value_or( 2500000 );
        min_len    = p.get_optional< std::size_t >( "min_len"    ).value_or( 0 );
        max_len    = p.get_optional< std::size_t >( "max_len"    ).value_or( 0 );
        node_path  = p.get_optional< std::string >( "node_path"  ).value_or( "/home/joyel/bin/node" );
        heatbub_js = p.get_optional< std::string >( "heatbub_js" ).value_or( "/home/joyel/WorkDir/AgoD3/heatmap_bubble_plot/heatmap_bubble_plot.js" );

        if( p.get_child_optional( "biotype_list"  )) emplace_list( p, "biotype_list", biotype_list );
        if( p.get_child_optional( "analysis_list" )) emplace_list( p, "analysis_list", analysis_list );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        auto& bed_samples  = db.bed_samples;
        auto& genome_table = db.genome_table;

        std::size_t analysis_size = analysis_list.empty() ? biotype_list.size() : analysis_list.size() ;

        monitor.set_monitor( "Component GeneTypeAnalyzer", 4 + analysis_size + ( output_annobed ? bed_samples.size() : 0 ));
        monitor.log( "Component GeneTypeAnalyzer", "Start" );


        monitor.log( "Component GeneTypeAnalyzer", "Filtering ... " );
        filtering( bed_samples, biotype_list, is_keep_other_biotype, max_anno_merge_size );
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        std::string output_path = db.output_dir().string() + ( db.output_dir().string().at( db.output_dir().string().length() -1 ) != '/' ? "/" : "" ) ;
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Biotypes" ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other_Seed" ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other" ));

        monitor.log( "Component GeneTypeAnalyzer", "Outputing ... Biotypes" );
        algorithm::GeneTypeAnalyzerBiotype( output_path + "Biotypes/", genome_table, bed_samples, biotype_list, min_len, max_len, sudo_count, webpage_update_only );
        algorithm::AnnoLengthIndexType ano_len_idx;

        std::vector< std::vector< algorithm::CountingTableType >> anno_table_tail;
        std::vector< std::map< std::string, std::string >> anno_mark;
        std::vector< std::map< std::string, std::map< std::string, double >>> seed_match_table;

        for( std::size_t i = 0; i < analysis_size; ++i )
        {
            // break;

            auto& biotype = analysis_list.empty() ? biotype_list[i] : analysis_list[i];
            if( biotype == "rmsk" ) continue;

            monitor.set_monitor( "\tBiotype Analysis - " + biotype, 6 );

            monitor.log( "Component GeneTypeAnalyzer", "Outputing ... " + biotype + " [ " + std::to_string( i+1 ) + " / " + std::to_string( analysis_size ) + " ]" );
            monitor.log( "\tBiotype Analysis - " + biotype, "Start" );

            anno_table_tail = std::vector< std::vector< algorithm::CountingTableType >>(
                    bed_samples.size(), std::vector< algorithm::CountingTableType >( 6, algorithm::CountingTableType() ));

            seed_match_table = std::vector< std::map< std::string, std::map< std::string, double >>>( bed_samples.size(), std::map< std::string, std::map< std::string, double >>());
            anno_mark = std::vector< std::map< std::string, std::string >>( bed_samples.size(), std::map< std::string, std::string >());

            monitor.log( "\tBiotype Analysis - " + biotype, "Makeing Annotation Counting Table ..." );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                smp_parallel_pool.job_post([ smp, &bed_samples, &anno_table_tail, &anno_mark, &genome_table, &biotype, this ] ()
                {
                    make_anno_table( bed_samples[ smp ].second, anno_table_tail[ smp ], anno_mark[ smp ], genome_table, biotype );
                });
            }

            smp_parallel_pool.flush_pool();

            ano_len_idx = get_ano_len_idx( genome_table, bed_samples, biotype );
            table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );

            std::size_t total_anno_counts = ano_len_idx.first.size();

            if( total_anno_counts == 0 )
            {
                monitor.log( "\tBiotype Analysis - " + biotype, "Skip with " + std::to_string( total_anno_counts ) + " annotations" );
                monitor.log( "\tBiotype Analysis - " + biotype, "Skip with " + std::to_string( total_anno_counts ) + " annotations" );
                monitor.log( "\tBiotype Analysis - " + biotype, "Skip with " + std::to_string( total_anno_counts ) + " annotations" );
                monitor.log( "\tBiotype Analysis - " + biotype, "Skip with " + std::to_string( total_anno_counts ) + " annotations" );
                continue;
            }

            monitor.log( "\tBiotype Analysis - " + biotype, "Doing Annotation Analysis ..." );

            if( !webpage_update_only )
            {
                if( biotype != "miRNA_mirtron" )
                {
                    boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other/" + biotype ));
                    boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other_Seed/" + biotype ));
                }
                else
                {
                    boost::filesystem::create_directory( boost::filesystem::path( output_path + "miR" ));
                    boost::filesystem::create_directory( boost::filesystem::path( output_path + "miR_Seed" ));
                }

                // if( biotype == "miRNA_mirtron" )
                {
                    // algorithm::GeneTypeAnalyzerQuantile( ano_len_idx, anno_table_tail );

                    std::vector< double > tmm_means;
                    algorithm::GeneTypeAnalyzerTmm( ano_len_idx, anno_table_tail, tmm_means, 30, 5, 2000 ); // Mg / Ag / mean

                    std::cerr << biotype;
                    for( auto& tmm : tmm_means ) std::cerr << "\t" << tmm;
                    std::cerr << "\n";
                }
            }

            algorithm::GeneTypeAnalyzerEachtype(
                      biotype
                    , output_path + ( biotype == "miRNA_mirtron" ? "miR/" : "Other/" )
                    , bed_samples
                    , ano_len_idx
                    , anno_table_tail
                    , seed_match_table
                    , anno_mark
                    , genome_table
                    , thread_number
                    , node_path
                    , heatbub_js
                    , min_len
                    , max_len
                    , extend_merge
                    , extend_refseq
                    , max_anno_merge_size
                    , webpage_update_only
                    , false
                    );

            monitor.log( "\tBiotype Analysis - " + biotype, "Makeing Seed Counting Table ..." );

            anno_table_tail = std::vector< std::vector< algorithm::CountingTableType >>(
                    bed_samples.size(), std::vector< algorithm::CountingTableType >( 6, algorithm::CountingTableType() ));

            anno_mark = std::vector< std::map< std::string, std::string >>( bed_samples.size(), std::map< std::string, std::string >());

            if( !webpage_update_only )
            {
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                {
                    smp_parallel_pool.job_post([ smp, &bed_samples, &anno_table_tail, &seed_match_table, &genome_table, &biotype, this ] ()
                    {
                        make_seed_table( bed_samples[ smp ].second, anno_table_tail[ smp ], seed_match_table[ smp ], genome_table, biotype );
                    });
                }

                smp_parallel_pool.flush_pool();

                ano_len_idx = get_ano_len_idx( genome_table, bed_samples, biotype, true );
                table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );
                seed_refinding( ano_len_idx, anno_table_tail, seed_match_table );

                // if( biotype == "miRNA_mirtron" )
                    algorithm::GeneTypeAnalyzerQuantile( ano_len_idx, anno_table_tail );
            }

            monitor.log( "\tBiotype Analysis - " + biotype, "Doing Seed Analysis ..." );

            algorithm::GeneTypeAnalyzerEachtype(
                      biotype
                    , output_path + ( biotype == "miRNA_mirtron" ? "miR_Seed/" : "Other_Seed/" )
                    , bed_samples
                    , ano_len_idx
                    , anno_table_tail
                    , seed_match_table
                    , anno_mark
                    , genome_table
                    , thread_number
                    , node_path
                    , heatbub_js
                    , min_len
                    , max_len
                    , extend_merge
                    , extend_refseq
                    , max_anno_merge_size
                    , webpage_update_only
                    , true
                    );

            monitor.log( "\tBiotype Analysis - " + biotype, "Complete with " + std::to_string( total_anno_counts ) + " annotations" );
        }

        if( output_annobed )
        {
            for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                monitor.log( "Component GeneTypeAnalyzer", "Outputing ... AnnoBed [ " + std::to_string( smp+1 ) + " / " + std::to_string( bed_samples.size() ) + " ]" );

                if( !webpage_update_only )
                {
                    std::ofstream annobed_output( output_path + bed_samples[ smp ].first + "_analyzedbed.text" );
                    annobed_outputing( annobed_output, genome_table, bed_samples[ smp ].second );
                    annobed_output.close();
                }
            }
        }

        monitor.log( "Component GeneTypeAnalyzer", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
