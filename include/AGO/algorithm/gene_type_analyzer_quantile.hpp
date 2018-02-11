#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/quantile.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerQuantile
{
    ParaThreadPool smp_parallel_pool;
    std::vector< QuantileDataType<> > ppm_qvec;
    std::vector< std::vector< CountingTableType >> rcvr_table_tail;

  public:

    GeneTypeAnalyzerQuantile()
        : smp_parallel_pool( 0 )
        , ppm_qvec( 0, QuantileDataType<>() )
        , rcvr_table_tail( 0, std::vector< CountingTableType >( 6, CountingTableType() ))
    {}

    GeneTypeAnalyzerQuantile(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::vector< CountingTableType >>& anno_table_tail
            )
        : smp_parallel_pool( anno_table_tail.size() )
        , ppm_qvec( anno_table_tail.size(), QuantileDataType<>() )
        , rcvr_table_tail( anno_table_tail.size(), std::vector< CountingTableType >( 6, CountingTableType() ))
    {
        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &ano_len_idx, &anno_table_tail, this ] ()
            {
                double sum = 0.0;
                std::vector< double > ppm_vec;

                for( auto& anno : ano_len_idx.first )
                {
                    sum = get_sum( anno, anno_table_tail[ smp ] );
                    make_recover(  anno, anno_table_tail[ smp ], rcvr_table_tail[ smp ], sum );
                    ppm_vec.emplace_back( sum );
                }

                ppm_qvec[ smp ] = QuantileDataType<>( ppm_vec );
            });
        }

        smp_parallel_pool.flush_pool();
        QuantileNor quntile( ppm_qvec );

        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &ano_len_idx, &anno_table_tail, this ] ()
            {
                recovering( ano_len_idx, anno_table_tail[ smp ], rcvr_table_tail[ smp ], ppm_qvec[ smp ] );
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
            QuantileDataType<>& ppm_qtl
            )
    {

        std::size_t idx = 0;

        for( auto& anno : ano_len_idx.first )
        {
            for( std::size_t tail = 0; tail < anno_table_tail.size(); ++tail )
                for( auto& len : anno_table_tail[ tail ][ anno ] )
                    len.second = rcvr_table_tail[ tail ][ anno ][ len.first ] * ppm_qtl.value_[ idx ].first;
            idx++;
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
