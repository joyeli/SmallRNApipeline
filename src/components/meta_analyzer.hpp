#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/meta_analyzer.hpp>

namespace ago {
namespace component {

class MetaAnalyzer : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component MetaAnalyzer", 1 );

        std::string output_path( db.output_dir().string() );

        std::vector< std::string > idxlen( get_index( db.analyzer_result_samples, ".LenDist_GMPM_ppm" ));
        std::vector< std::string > idxmir( get_index( db.analyzer_result_samples, ".MirDist_GMPM_ppm" ));

        std::vector< QuantileDataType<> > lenq = quantile( db.analyzer_result_samples, idxlen, ".LenDist_GMPM_ppm" );
        std::vector< QuantileDataType<> > mirq = quantile( db.analyzer_result_samples, idxmir, ".MirDist_GMPM_ppm" );

        ago::algorithm::MetaAnalyzer meta;

        monitor.set_monitor( "MetaAnalyzing", 10 );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , lenq
            , idxlen
            , ".LenDist"
            , "GMPM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , lenq
            , idxlen
            , ".LenDist"
            , "GM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , lenq
            , idxlen
            , ".LenDist"
            , "PM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.tailing_ratio(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , ".LenDist_Tailing_Ratio"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , mirq
            , idxmir
            , ".MirDist"
            , "GMPM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , mirq
            , idxmir
            , ".MirDist"
            , "GM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.quantile_transfer(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , mirq
            , idxmir
            , ".MirDist"
            , "PM"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.tailing_ratio(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , ".MirDist_Tailing_Ratio"
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.tailing_ratio_for_pm(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , ".PM_Tailing_Ratio"
            , db.quantile_result_samples[ 2 ].second   // lendist_gm_readcount
            , db.quantile_result_samples[ 3 ].second   // lendist_gm_ppm
            , db.quantile_result_samples[ 4 ].second   // lendist_pm_readcount
            , db.quantile_result_samples[ 5 ].second   // lendist_pm_ppm
        );

        monitor.log( "MetaAnalyzing", " ... " );

        meta.tailing_ratio_for_mir(
              db.quantile_result_samples
            , db.analyzer_result_samples
            , ".miRNA_Tailing_Ratio"
            , db.quantile_result_samples[ 8 ].second   // mirdist_gmpm_ppm
        );

        monitor.log( "MetaAnalyzing", " ... Complete" );

        db.analyzer_result_samples.clear();

        monitor.log( "Component MetaAnalyzer", "Complete!!" );
    }

    std::vector< std::string > get_index(
        std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
        , const char* ppm
    )
    {
        std::vector< std::string > idx;

        for( auto& sample : analyzer_result_samples )
        {
            for( auto& tables : sample.second )
            {
                if( tables.find( ppm ) != tables.end() )
                {
                    for( auto& pair : tables )
                    {
                        if( pair.first == ".LenDist_GMPM_ppm" ||
                                pair.first == ".MirDist_GMPM_ppm" ||
                                pair.first == "SUM_ANNO_SUM_ANNO" ||
                                pair.first == "SUM_ANNO" )
                        {}
                        else
                        {
                            idx.emplace_back( pair.first );
                        }
                    }
                }
            }
        }

        std::sort( idx.begin(), idx.end() );
        auto it = std::unique( idx.begin(), idx.end() );
        idx.resize( std::distance( idx.begin(), it ));

        return idx;
    }

    std::vector< QuantileDataType<> > quantile(
          std::vector< std::pair< std::string, ago::engine::DataPool::AnalyzerResultType >>& analyzer_result_samples
        , const std::vector< std::string >& index
        , const char* ppm
    )
    {
        std::vector< QuantileDataType<> > ppms_qvec;

        for( auto& sample : analyzer_result_samples )
        {
            for( auto& annos : sample.second )
            {
                if( annos.find( ppm ) != annos.end() )
                {
                    std::vector< double > ppms;

                    for( auto& idx : index )
                    {
                        auto it = annos.find( idx );

                        if( it != annos.end() )
                        {
                            ppms.push_back( it->second[ "SUM_LEN" ] );
                        }
                        else
                        {
                            ppms.push_back( 0 );
                        }
                    }

                    ppms_qvec.push_back( QuantileDataType<>( ppms ));
                }
            }
        }

        QuantileNor quntile( ppms_qvec );
        return ppms_qvec;
    }

};

} // end of namespace component
} // end of namespace ago
