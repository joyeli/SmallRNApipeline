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
    std::string rnafold_path;

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
        rnafold_path = p.get_optional< std::string >( "rnafold_path"  ).value_or( "" );

        if( p.get_child_optional( "biotype_list"  )) emplace_list( p, "biotype_list", biotype_list );
        if( p.get_child_optional( "analysis_list" )) emplace_list( p, "analysis_list", analysis_list );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {}

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        auto& bed_samples  = db.bed_samples;
        auto& genome_table = db.genome_table;

        std::vector< double > tmm_means;
        std::vector< std::pair< std::string, std::vector< double >>> trim_means;

        std::size_t analysis_size = analysis_list.empty() ? biotype_list.size() : analysis_list.size() ;

        monitor.set_monitor( "Component GeneTypeAnalyzer", 4 + analysis_size );
        monitor.log( "Component GeneTypeAnalyzer", "Start" );

        monitor.log( "Component GeneTypeAnalyzer", "Filtering ... " );
        filtering( bed_samples, biotype_list, is_keep_other_biotype, max_anno_merge_size );

        std::string output_path = db.output_dir().string() + ( db.output_dir().string().at( db.output_dir().string().length() -1 ) != '/' ? "/" : "" ) ;

        monitor.log( "Component GeneTypeAnalyzer", "Outputing ... Biotypes" );

        boost::filesystem::create_directory( boost::filesystem::path( output_path + "Biotypes" ));
        algorithm::GeneTypeAnalyzerBiotype( output_path + "Biotypes/", genome_table, bed_samples, biotype_list, min_len, max_len, sudo_count, webpage_update_only );

        algorithm::AnnoLengthIndexType ano_len_idx;

        std::vector< std::vector< algorithm::CountingTableType >> anno_table_tail;
        std::vector< std::map< std::string, std::string >> anno_mark;

        ParaThreadPool ana_parallel_pool( thread_number );

        for( std::size_t i = 0; i < analysis_size; ++i )
        {
            // break;

            auto& biotype = analysis_list.empty() ? biotype_list[i] : analysis_list[i];
            if( biotype == "rmsk" ) continue;

            monitor.set_monitor( "\tBiotype Analysis - " + biotype, 2 );
            monitor.log( "Component GeneTypeAnalyzer", "Outputing ... " + biotype + " [ " + std::to_string( i+1 ) + " / " + std::to_string( analysis_size ) + " ]" );
            monitor.log( "\tBiotype Analysis - " + biotype, "Start" );

            anno_table_tail = std::vector< std::vector< algorithm::CountingTableType >>( bed_samples.size(), std::vector< algorithm::CountingTableType >( 6 ));
            anno_mark = std::vector< std::map< std::string, std::string >>( bed_samples.size() );
            ano_len_idx = get_ano_len_idx( genome_table, bed_samples, biotype );

            if( ano_len_idx.first.empty() )
            {
                monitor.log( "\tBiotype Analysis - " + biotype, "Skip with 0 annotations" );
                continue;
            }

            if( !webpage_update_only )
            {
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                {
                    ana_parallel_pool.job_post([ smp, &bed_samples, &anno_table_tail, &anno_mark, &genome_table, &biotype, this ] ()
                    {
                        make_anno_table( bed_samples[ smp ].second, anno_table_tail[ smp ], anno_mark[ smp ], genome_table, biotype );
                    });
                }

                ana_parallel_pool.flush_pool();

                table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );
                algorithm::GeneTypeAnalyzerQuantile( ano_len_idx, anno_table_tail );

                // algorithm::GeneTypeAnalyzerTmm( ano_len_idx, anno_table_tail, tmm_means, 30, 5, 2000 ); // Mg / Ag / mean
                // table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );
                // trim_means.emplace_back( std::make_pair( biotype, tmm_means ));
                // tmm_means.clear();
            }

            ana_parallel_pool.job_post([ &biotype, &bed_samples, &genome_table, &anno_table_tail, &ano_len_idx, &anno_mark, &output_path, this ] ()
            {
                do_analysis( biotype, bed_samples, genome_table, anno_table_tail, ano_len_idx, anno_mark, output_path );
            });
            ana_parallel_pool.job_post([ &biotype, &bed_samples, &genome_table, &anno_table_tail, &ano_len_idx, &anno_mark, &output_path, this ] ()
            {
                do_analysis( biotype, bed_samples, genome_table, anno_table_tail, ano_len_idx, anno_mark, output_path, "Seed" );
            });

            for( std::size_t len = min_len; len <= max_len; ++len )
                ana_parallel_pool.job_post([ len, &biotype, &bed_samples, &genome_table, &anno_table_tail, &ano_len_idx, &anno_mark, &output_path, this ] ()
                {
                    do_analysis( biotype, bed_samples, genome_table, anno_table_tail, ano_len_idx, anno_mark, output_path, std::to_string( len ));
                });

            if( biotype == "miRNA_mirtron" || biotype == "miRNA" || biotype == "mirtron" )
            {
                ana_parallel_pool.job_post([ &biotype, &bed_samples, &genome_table, &anno_table_tail, &ano_len_idx, &anno_mark, &output_path, this ] ()
                {
                    do_analysis( biotype, bed_samples, genome_table, anno_table_tail, ano_len_idx, anno_mark, output_path, "5p" );
                });
                ana_parallel_pool.job_post([ &biotype, &bed_samples, &genome_table, &anno_table_tail, &ano_len_idx, &anno_mark, &output_path, this ] ()
                {
                    do_analysis( biotype, bed_samples, genome_table, anno_table_tail, ano_len_idx, anno_mark, output_path, "3p" );
                });
            }

            ana_parallel_pool.flush_pool();

            monitor.log( "\tBiotype Analysis - " + biotype, "Complete with " + std::to_string( ano_len_idx.first.size() ) + " annotations" );
        }

        if( !webpage_update_only )
        {
            if( output_annobed )
            {
                for( size_t smp = 0; smp < bed_samples.size(); ++smp )
                {
                    std::ofstream annobed_output( output_path + bed_samples[ smp ].first + "_analyzedbed.text" );
                    annobed_outputing( annobed_output, genome_table, bed_samples[ smp ].second );
                    annobed_output.close();
                }
            }

            // std::ofstream tmm_output( output_path + "tmm_means.text" );
            // tmm_output << "Biotypes";

            // for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            //     tmm_output << "\t" << bed_samples[ smp ].first;

            // for( auto& bio : trim_means )
            // {
            //     tmm_output << "\n" << bio.first;

            //     for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            //         tmm_output << "\t" << bio.second[ smp ];
            // }

            // tmm_output << "\n";
            // tmm_output.close();
        }

        monitor.log( "Component GeneTypeAnalyzer", "Complete" );
    }

    void do_analysis( auto& biotype, auto& bed_samples, auto& genome_table, auto anno_table_tail, auto ano_len_idx, auto anno_mark, auto output_path, const std::string& token = "" )
    {
        output_path = output_path + "Genetype" + ( token == "" ? "" : ( "_" + token )) + "/";

        boost::filesystem::create_directory( boost::filesystem::path( output_path ));
        boost::filesystem::create_directory( boost::filesystem::path( output_path + biotype ));

        std::vector< std::map< std::string, std::map< std::string, double >>> seed_match_table
            = std::vector< std::map< std::string, std::map< std::string, double >>>( bed_samples.size() );

        if( token != "" && !webpage_update_only )
        {
            ParaThreadPool smp_parallel_pool( bed_samples.size() +1 );

            smp_parallel_pool.job_post([ &ano_len_idx, &genome_table, &bed_samples, &biotype, &token, this ] ()
            {
                ano_len_idx = get_ano_len_idx( genome_table, bed_samples, biotype, token );
            });

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                smp_parallel_pool.job_post([ smp, &anno_table_tail, &seed_match_table, &anno_mark, &genome_table, &biotype, &token, this ] ()
                {
                    make_ana_table( anno_table_tail[ smp ], anno_mark[ smp ], seed_match_table[ smp ], token );
                });
            }

            smp_parallel_pool.flush_pool();

            table_refinding(
                  ano_len_idx
                , anno_table_tail
                , ( token != "Seed" & token != "" ? std::stoi( token ) : min_len )
                , ( token != "Seed" & token != "" ? std::stoi( token ) : max_len )
                , sudo_count
                );

            if( token == "Seed" )
                seed_refinding( ano_len_idx, anno_table_tail, seed_match_table );
        }

        algorithm::GeneTypeAnalyzerEachtype(
                  biotype
                , output_path
                , bed_samples 
                , ano_len_idx
                , anno_table_tail
                , seed_match_table
                , anno_mark
                , genome_table
                , thread_number
                , node_path
                , heatbub_js
                , ( token != "Seed" & token != "" ? std::stoi( token ) : min_len )
                , ( token != "Seed" & token != "" ? std::stoi( token ) : max_len )
                , extend_merge
                , extend_refseq
                , max_anno_merge_size
                , webpage_update_only
                , rnafold_path
                , token
                );
    }
};

} // end of namespace component
} // end of namespace ago
