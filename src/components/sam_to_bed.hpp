#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/converter/sam2rawbed.hpp>
#include <mutex>

namespace ago {
namespace component {

class SamToBed : public engine::NamedComponent
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

        monitor.set_monitor( "Component SamToBed", db.sam_samples.size() +2 );
        monitor.log( "Component SamToBed", "Start" );

        for( auto& sam_sample : db.sam_samples )
        {
            for( auto sam = sam_sample.second.begin(); sam != sam_sample.second.end(); ++sam )
            {
                std::get<3>( sam->data ) = std::get<3>( sam->data ) - 1;
            }

            Sam2RawBed< std::vector< Sam<> >* > sam2bed;
            auto raw_beds( sam2bed.Convert( &sam_sample.second ));

            std::vector< AnnotationRawBed<> > annotation_rawbeds;

            for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
            {
                annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
            }

            sam2bed.rawbed_map_->clear();
            sam2bed.rawbed_map2_->clear();

            db.bed_samples.emplace_back( sam_sample.first, annotation_rawbeds );
            monitor.log( "Component SamToBed", ( sam_sample.first ).c_str() );
        }

        db.sam_samples.clear();
        monitor.log( "Component SamToBed", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
