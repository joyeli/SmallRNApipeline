#pragma once
#include <AGO/engine/components/named_component.hpp>
// #include <pokemon/converter/sam2rawbed.hpp>
#include <AGO/format/md_rawbed.hpp>
#include <mutex>

namespace ago {
namespace component {

class SamToBed : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    std::size_t max_tail_len;

  public:

    using Base::Base;

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        max_tail_len = p.get_optional< std::size_t >( "max_tail_len" ).value_or( 5 );
    }

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component SamToBed", db.sam_samples.size() +2 );
        monitor.log( "Component SamToBed", "Start" );

        std::map< ago::format::MDRawBed, std::size_t > md_rawbeds_map;
        std::vector< ago::format::MDRawBed > md_rawbeds;
        ago::format::MDRawBed md_rawbed;

        for( auto& sam_sample : db.sam_samples )
        {
            // Sam2RawBed< std::vector< Sam<> >* > sam2bed;
            // auto raw_beds( sam2bed.Convert( &sam_sample.second ));

            // std::vector< AnnotationRawBed<> > annotation_rawbeds;

            for( auto& sam : sam_sample.second )
            {
                md_rawbed = ago::format::MDRawBed( sam );
                md_rawbed.reducing_tail( max_tail_len );

                if( md_rawbeds_map.find( md_rawbed ) == md_rawbeds_map.end() )
                    md_rawbeds_map[ md_rawbed ] = 0;

                md_rawbeds_map[ md_rawbed ]++;
            }

            for( auto& mdbed : md_rawbeds_map )
            {
                (*(( uint32_t* )( mdbed.first.get_reads_count() ))) = mdbed.second;
                md_rawbeds.emplace_back( mdbed.first );
            }

            // for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
            // {
            //     annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
            // }

            // sam2bed.rawbed_map_->clear();
            // sam2bed.rawbed_map2_->clear();

            // db.bed_samples.emplace_back( sam_sample.first, annotation_rawbeds );

            md_rawbeds_map.clear();
            db.bed_samples.emplace_back( sam_sample.first, md_rawbeds );
            monitor.log( "Component SamToBed", ( sam_sample.first ).c_str() );
            md_rawbeds.clear();
        }

        db.sam_samples.clear();
        monitor.log( "Component SamToBed", "Complete" );
    }
};

} // end of namespace component
} // end of namespace ago
