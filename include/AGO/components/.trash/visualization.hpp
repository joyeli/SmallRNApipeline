#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/visualization.hpp>

namespace ago {
namespace component {

class Visualization : public engine::NamedComponent
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

        monitor.set_monitor( "Component Visualization", 9 );
        monitor.log( "Component Visualization", "Start" );

        std::string output_path = db.output_dir().string();
        ago::algorithm::Visualization visualization;

		visualization.make_html( output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_biotype( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_lendist( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

	  	visualization.make_mirdist( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_dotplot( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_valplot( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_lenplus( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", " ... " );

		visualization.make_mirtail( db.quantile_result_samples, output_path );

        monitor.log( "Component Visualization", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
