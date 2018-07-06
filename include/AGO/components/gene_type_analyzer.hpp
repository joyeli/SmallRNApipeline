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
{
    using Base = engine::NamedComponent;

    std::vector< std::string > biotype_list;

    double sudo_count;

    bool output_annobed;
    bool is_keep_other_biotype;

    std::size_t thread_number;
    std::size_t extend_refseq;
    std::size_t extend_merge;
    std::size_t min_len;
    std::size_t max_len;

    std::string node_path;
    std::string heatbub_js;

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        sudo_count = p.get_optional< double   >( "sudo_count"     ).value_or( 0.000001 );
        output_annobed = p.get_optional< bool >( "output_annobed" ).value_or( true     );
        is_keep_other_biotype = p.get_optional< bool >( "is_keep_other_biotype" ).value_or( false );
        thread_number  = p.get_optional< std::size_t >( "thread_number" ).value_or( 8  );
        extend_refseq  = p.get_optional< std::size_t >( "extend_refseq" ).value_or( 10 );
        extend_merge   = p.get_optional< std::size_t >( "extend_merge"  ).value_or( 2  );
        min_len    = p.get_optional< std::size_t >( "min_len"    ).value_or( 0 );
        max_len    = p.get_optional< std::size_t >( "max_len"    ).value_or( 0 );
        node_path  = p.get_optional< std::string >( "node_path"  ).value_or( "/home/joyel/bin/node" );
        heatbub_js = p.get_optional< std::string >( "heatbub_js" ).value_or( "/home/joyel/WorkDir/AgoD3/heatmap_bubble_plot/heatmap_bubble_plot.js" );

        if(  p.get_child_optional( "biotype_list" ))
        {
            std::size_t i = 0;
            std::size_t miRNA_idx = 0;
            std::size_t mirtron_idx = 0;
            std::size_t miRNA_mirtron_idx = 0;

            for( auto& biotype : p.get_child( "biotype_list" ))
            {
                i++;
                if( biotype.second.data() == "miRNA" ) miRNA_idx = i; 
                if( biotype.second.data() == "mirtron" ) mirtron_idx = i; 
                if( biotype.second.data() == "miRNA_mirtron" ) miRNA_mirtron_idx = i; 
                biotype_list.emplace_back( biotype.second.data() );
            }

            if(( miRNA_mirtron_idx != 0 && ( miRNA_idx == 0 || mirtron_idx == 0 ))
            || ( miRNA_idx != 0 && mirtron_idx != 0 && miRNA_idx < mirtron_idx  ))
            {
                std::vector< std::string > biotype_list_temp;

                for( std::size_t i = 0; i < biotype_list.size(); ++i )
                {
                    if( biotype_list[i] == "mirtron" ) continue;
                    if( miRNA_mirtron_idx != 0 )
                    {
                        if( biotype_list[i] == "miRNA" ) continue;
                        if( biotype_list[i] == "miRNA_mirtron" )
                        {
                            biotype_list_temp.emplace_back( "mirtron" );
                            biotype_list_temp.emplace_back( "miRNA" );
                        }
                    }
                    else if( biotype_list[i] == "miRNA" )
                    {
                        biotype_list_temp.emplace_back( "mirtron" );
                        biotype_list_temp.emplace_back( "miRNA" );
                    }

                    biotype_list_temp.emplace_back( biotype_list[i] );
                }

                biotype_list = biotype_list_temp;
            }

            biotype_list.emplace_back( "un_annotated" );
        }
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

        monitor.set_monitor( "Component GeneTypeAnalyzer", 5 + biotype_list.size() + ( output_annobed ? bed_samples.size() : 0 ));
        monitor.log( "Component GeneTypeAnalyzer", "Start" );


        monitor.log( "Component GeneTypeAnalyzer", "Filtering ... " );
        filtering( bed_samples, biotype_list, is_keep_other_biotype );
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        std::string output_path = db.output_dir().string() + ( db.output_dir().string().at( db.output_dir().string().length() -1 ) != '/' ? "/" : "" ) ;
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Biotypes" ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other" ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "miR" ));

        monitor.log( "Component GeneTypeAnalyzer", "Outputing ... Biotypes" );
        algorithm::GeneTypeAnalyzerBiotype( output_path + "Biotypes/", genome_table, bed_samples, biotype_list, min_len, max_len, sudo_count );
        algorithm::AnnoLengthIndexType ano_len_idx;

        monitor.log( "Component GeneTypeAnalyzer", "Making ... Counting Table" );
        std::vector< std::vector< algorithm::CountingTableType >> anno_table_tail;
        std::vector< std::map< std::string, std::string >> anno_mark;

        for( std::size_t i = 0; i < biotype_list.size(); ++i )
        {
            auto& biotype = biotype_list[i];
            monitor.log( "Component GeneTypeAnalyzer", "Outputing ... " + biotype + " [ " + std::to_string( i+1 ) + " / " + std::to_string( biotype_list.size() ) + " ]" );

            if( biotype == "rmsk" )
                continue;

            if( biotype != "miRNA_mirtron" )
                boost::filesystem::create_directory( boost::filesystem::path( output_path + "Other/" + biotype ));

            anno_table_tail = std::vector< std::vector< algorithm::CountingTableType >>(
                    bed_samples.size(), std::vector< algorithm::CountingTableType >( 6, algorithm::CountingTableType() ));

            anno_mark = std::vector< std::map< std::string, std::string >>( bed_samples.size(), std::map< std::string, std::string >());

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

            if( biotype == "miRNA_mirtron" )
                algorithm::GeneTypeAnalyzerQuantile( ano_len_idx, anno_table_tail );

            algorithm::GeneTypeAnalyzerEachtype(
                      biotype
                    , output_path + ( biotype == "miRNA_mirtron" ? "miR/" : "Other/" )
                    , bed_samples
                    , ano_len_idx
                    , anno_table_tail
                    , anno_mark
                    , thread_number
                    , genome_table
                    , node_path
                    , heatbub_js
                    , min_len
                    , max_len
                    , extend_merge
                    , extend_refseq
                    );
        }

        if( output_annobed )
        {
            for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                monitor.log( "Component GeneTypeAnalyzer", "Outputing ... AnnoBed [ " + std::to_string( smp+1 ) + " / " + std::to_string( bed_samples.size() ) + " ]" );

                std::ofstream annobed_output( output_path + bed_samples[ smp ].first + "_analyzedbed.text" );
                annobed_outputing( annobed_output, genome_table, bed_samples[ smp ].second );
                annobed_output.close();
            }
        }

        monitor.log( "Component GeneTypeAnalyzer", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
