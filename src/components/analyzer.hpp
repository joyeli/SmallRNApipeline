#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/analyzer_lendist.hpp>
#include <AGO/algorithm/analyzer_mirdist.hpp>

namespace ago {
namespace component {

class Analyzer : public engine::NamedComponent
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

        monitor.set_monitor( "Component Analyzer", db.bed_samples.size() +2);
        monitor.log( "Component Analyzer", "Start" );

        std::string output_path( db.output_dir().string() );

        ago::algorithm::AnalyzerLenDist len;
        ago::algorithm::AnalyzerMirDist mir;

        for( auto& sample : db.bed_samples )
        {

            len.analyzer.run( &sample.second, 0, true, output_path, sample.first, db.genome_table, db.analyzer_result, 0 );
            len.tailing_ratio( db.analyzer_result );

            mir.analyzer.run( &sample.second, 0, true, output_path, sample.first, db.genome_table, db.analyzer_result, 0 );
            mir.tailing_ratio( db.analyzer_result );
            mir.tailing_ratio_for_PM( db.analyzer_result, sample.second, db.genome_table );
            mir.tailing_ratio_for_miRNA( db.analyzer_result, sample.second, db.genome_table );
            mir.tailing_detail_for_miRNA( sample, db.genome_table, output_path );

            monitor.log( "Component Analyzer", ( sample.first ).c_str() );

            db.analyzer_result_samples.emplace_back( sample.first, db.analyzer_result );
            db.analyzer_result.clear();
            sample.second.clear();
        }

        db.bed_samples.clear();

        monitor.log( "Component Analyzer", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
