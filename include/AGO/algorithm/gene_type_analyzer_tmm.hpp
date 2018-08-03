#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>
#include <AGO/algorithm/tmm_for_pipeline.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerTmm
{
    ParaThreadPool smp_parallel_pool;
    std::vector< TmmDataType<> > tmm_vec;
    std::vector< std::vector< CountingTableType >> rcvr_table_tail;
    std::vector< std::map< std::string, std::map< std::string, double >>> rcvr_table_anno;
    std::vector< std::map< std::string, std::size_t >> anno_idx;

  public:

    GeneTypeAnalyzerTmm()
        : smp_parallel_pool( 0 )
        , tmm_vec( 0 )
        , rcvr_table_tail( 0 )
        , rcvr_table_anno( 0 )
        , anno_idx( 0 )
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
        , tmm_vec( anno_table_tail.size() )
        , rcvr_table_tail( anno_table_tail.size(), std::vector< CountingTableType >( 6, CountingTableType() ))
        , rcvr_table_anno( anno_table_tail.size() )
        , anno_idx( anno_table_tail.size() )
    {
        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &ano_len_idx, &anno_table_tail, this ] ()
            {
                double sum = 0.0;
                std::size_t idx = 0;
                std::string anno_name;

                std::vector< double > ppm_vec;
                std::map< std::string, double > ppm_map;

                for( auto& anno : ano_len_idx.first )
                {
                    sum = get_sum( anno, anno_table_tail[ smp ] );
                    make_recover(  anno, anno_table_tail[ smp ], rcvr_table_tail[ smp ], sum );

                    anno_name = get_anno_name( anno );
                    rcvr_table_anno[ smp ][ anno_name ][ anno ] = sum;

                    if( ppm_map.find( anno_name ) == ppm_map.end() )
                        ppm_map[ anno_name ] = 0.0;

                    ppm_map[ anno_name ] += sum;
                }

                for( auto& anno : ppm_map )
                {
                    anno_idx[ smp ][ anno.first ] = idx;
                    ppm_vec.emplace_back( anno.second );
                    idx++;

                    for( auto& anno_seed : rcvr_table_anno[ smp ][ anno.first ] )
                        anno_seed.second = anno_seed.second / anno.second;
                }

                tmm_vec[ smp ] = TmmDataType<>( ppm_vec );
            });
        }

        smp_parallel_pool.flush_pool();
        TmmNor TMM( tmm_vec, tmm_means, mg_per, ag_per, trim_mean );

        for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &ano_len_idx, &anno_table_tail, this ] ()
            {
                recovering( ano_len_idx, anno_table_tail[ smp ], rcvr_table_tail[ smp ], rcvr_table_anno[ smp ], anno_idx[ smp ], tmm_vec[ smp ] );
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
