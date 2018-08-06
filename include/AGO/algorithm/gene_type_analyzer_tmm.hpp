#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/tmm_for_pipeline.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerTmm
{
    ParaThreadPool smp_parallel_pool;

  public:

    GeneTypeAnalyzerTmm()
        : smp_parallel_pool( 0 )
    {}

    GeneTypeAnalyzerTmm(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail,
            std::vector< double >& tmm_means,
            const double& mg_per,
            const double& ag_per,
            const double& trim_mean
            )
        : smp_parallel_pool( anno_table_tail.size() )
    {
        tmm_means = std::vector< double >( anno_table_tail.size() );

        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &ano_len_idx, &anno_table_tail, &mg_per, &ag_per, &trim_mean, &tmm_means, this ] ()
            {
                double sum = 0.0;
                double tmm_mean = 0.0;

                std::size_t idx = 0;
                std::string anno_name;

                std::vector< double > ppm_vec;
                std::vector< TmmDataType<> > tmm_vec(1);

                std::map< std::string, double > ppm_map;
                std::map< std::string, std::size_t > anno_idx;

                std::vector< CountingTableType > rcvr_table_tail( 6 );
                std::map< std::string, std::map< std::string, double >> rcvr_table_anno;

                for( auto& anno : ano_len_idx.first )
                {
                    anno_name = get_anno_name( anno );
                    sum = get_sum( anno, anno_table_tail[ smp ] );

                    if( sum == 0 ) continue;

                    make_recover(  anno, anno_table_tail[ smp ], rcvr_table_tail, sum );
                    rcvr_table_anno[ anno_name ][ anno ] = sum;

                    if( ppm_map.find( anno_name ) == ppm_map.end() )
                        ppm_map[ anno_name ] = 0.0;

                    ppm_map[ anno_name ] += sum;
                }

                for( auto& anno : ppm_map )
                {
                    anno_idx[ anno.first ] = idx;
                    ppm_vec.emplace_back( anno.second );
                    idx++;

                    for( auto& anno_seed : rcvr_table_anno[ anno.first ] )
                        anno_seed.second = anno_seed.second / anno.second;
                }

                tmm_vec[0] = TmmDataType<>( ppm_vec );
                TmmNor TMM( tmm_vec, tmm_mean, mg_per, ag_per, trim_mean );

                tmm_means[ smp ] = tmm_mean;
                recovering( ano_len_idx, anno_table_tail[ smp ], rcvr_table_tail, rcvr_table_anno, anno_idx, tmm_vec[0] );
            });
        }

        smp_parallel_pool.flush_pool();
    }

    double get_sum( const std::string& anno, std::vector< CountingTableType >& anno_table_tail )
    {
        double res = 0.0;

        for( auto& anno_table : anno_table_tail )
            for( auto& len : anno_table[ anno ] )
                res += len.second;

        return res;
    }

    std::string get_anno_name( const std::string& anno )
    {
        std::vector< std::string > split;
        boost::iter_split( split, anno, boost::algorithm::first_finder( "_" ));

        std::string anno_name = split[0];

        for( std::size_t i = 1; i < split.size() -1; ++i )
            anno_name += "_" + split[i];

        return anno_name;
    }

    void make_recover(
            const std::string& anno,
            std::vector< CountingTableType >& anno_table_tail,
            std::vector< CountingTableType >& rcvr_table_tail,
            const double& sum
            )
    {
        for( std::size_t tail = 0; tail < anno_table_tail.size(); ++tail )
            for( auto& len : anno_table_tail[ tail ][ anno ] )
                rcvr_table_tail[ tail ][ anno ][ len.first ] = len.second / sum; 
    }
    
    void recovering(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< CountingTableType >& anno_table_tail,
            std::vector< CountingTableType >& rcvr_table_tail,
            std::map< std::string, std::map< std::string, double >>& rcvr_table_anno,
            std::map< std::string, std::size_t >& anno_idx,
            TmmDataType<>& tmm_ppms
            )
    {
        std::string anno_name;

        for( auto& anno : ano_len_idx.first )
        {
            anno_name = get_anno_name( anno );

            for( std::size_t tail = 0; tail < anno_table_tail.size(); ++tail )
                for( auto& len : anno_table_tail[ tail ][ anno ])
                    len.second
                        = rcvr_table_tail[ tail ][ anno ][ len.first ]
                        * rcvr_table_anno[ anno_name ][ anno ]
                        * tmm_ppms.values_[ anno_idx[ anno_name ]]
                        ;
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
