#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/gene_type_analyzer_filtering.hpp>
#include <AGO/algorithm/gene_type_analyzer_counting.hpp>
#include <AGO/algorithm/gene_type_analyzer_biotype.hpp>
#include <AGO/algorithm/gene_type_analyzer_quantile.hpp>
#include <AGO/algorithm/gene_type_analyzer_eachtype.hpp>
#include <AGO/algorithm/gene_type_analyzer_difference.hpp>
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
      , algorithm::GeneTypeAnalyzerDifference
      , algorithm::GeneTypeAnalyzerAnnobed
{
    using Base = engine::NamedComponent;

    std::vector< std::string > biotype_list;

    double sudo_count;

    bool is_filter_drop;
    bool output_annobed;

    std::size_t min_len;
    std::size_t max_len;

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        sudo_count     = p.get_optional< double >( "sudo_count"     ).value_or( 0.000001 );
        is_filter_drop = p.get_optional< bool   >( "is_filter_drop" ).value_or( true     );
        output_annobed = p.get_optional< bool   >( "output_annobed" ).value_or( true     );
        min_len = p.get_optional< std::size_t   >( "min_len"        ).value_or( 0        );
        max_len = p.get_optional< std::size_t   >( "max_len"        ).value_or( 0        );

        if(  p.get_child_optional( "biotype_list" ))
        {
            for( auto& biotype : p.get_child( "biotype_list" ))
            {
                biotype_list.emplace_back( biotype.second.data() );
            }
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

        monitor.set_monitor( "Component GeneTypeAnalyzer", 2 );
        monitor.log( "Component GeneTypeAnalyzer", "Start" );

        // auto bed_samples = db.bed_samples; for old style analysis
        auto& bed_samples  = db.bed_samples;
        auto& genome_table = db.genome_table;

        drop_filtering( bed_samples );
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        std::string output_path = db.output_dir().string() + ( db.output_dir().string().at( db.output_dir().string().length() -1 ) != '/' ? "/" : "" ) ;
        boost::filesystem::create_directory( boost::filesystem::path( output_path + "biotypes" ));

        algorithm::GeneTypeAnalyzerBiotype( output_path + "biotypes", genome_table, bed_samples, min_len, max_len, sudo_count );
        algorithm::AnnoLengthIndexType ano_len_idx;

        std::vector< std::vector< algorithm::CountingTableType >> anno_table_tail;
        std::vector< std::map< std::string, std::string >> anno_mark;

        for( std::size_t i = 0; i < biotype_list.size(); ++i )
        {
            auto& biotype = biotype_list[i];
            boost::filesystem::create_directory( boost::filesystem::path( output_path + biotype ));

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

            ano_len_idx = get_ano_len_idx( db.genome_table, bed_samples, biotype );
            table_refinding( ano_len_idx, anno_table_tail, min_len, max_len, sudo_count );

            algorithm::GeneTypeAnalyzerQuantile( ano_len_idx, anno_table_tail );
            algorithm::GeneTypeAnalyzerEachtype( output_path + biotype, bed_samples, ano_len_idx, anno_table_tail, anno_mark );
            algorithm::GeneTypeAnalyzerDifference( biotype, output_path, bed_samples, ano_len_idx, anno_table_tail );
        }

        if( output_annobed )
        {
            for( size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                std::ofstream annobed_output( output_path + bed_samples[ smp ].first + "_annobed.tsv" );
                annobed_outputing( annobed_output, db.genome_table, bed_samples[ smp ].second );
                annobed_output.close();
            }
        }

        monitor.log( "Component GeneTypeAnalyzer", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
