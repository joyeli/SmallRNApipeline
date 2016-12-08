#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class Annotator : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    using BedFileReaderImpl =  FileReader_impl<
          Bed
        , std::tuple< std::string, uint32_t, uint32_t, char, std::string, std::string >
        , SOURCE_TYPE::IFSTREAM_TYPE
    >;

    using AnnotationTrait = Annotation<
          BedFileReaderImpl
        , AnnoIgnoreStrand::NO_IGNORE
        , AnnoType::INTERSET
    >;

    using Annotations = AnnotationSet <
          AnnotationRawBed<>
        , AnnotationTrait
    >;

    std::vector< std::string > annotation_files_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "annotation_files" ))
        {
            db.push_path( "annotation_files", child.second );
        }
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );

        if( db.exist_path_tag( "annotation_files" ))
        {
            for( auto& bed : db.get_path_list( "annotation_files" ))
            {
                annotation_files_.emplace_back( bed.string() );
            }
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component Annotator", 1 );

        Annotations annotator( annotation_files_ );

        monitor.set_monitor( "Annotating", db.rawbed_samples.size()+1 );

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( db.rawbed_samples.size() );
        std::map< std::string, std::vector< AnnotationRawBed<> >> rawbed_samples_map;

        for( auto& sample : db.rawbed_samples )
        {
            smp_parallel_pool.job_post( [ sample, &annotator, &monitor, &smp_mutex, &rawbed_samples_map ] () mutable
            {
                for( auto& anno_rawbed : sample.second )
                {
                    annotator.AnnotateAll( anno_rawbed );
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );
                    rawbed_samples_map.emplace( sample ); 

                    monitor.log( "Annotating", " ... " + sample.first );
                }
            });
        }

        smp_parallel_pool.flush_pool();
        Annotations::clear_database();

        for( auto& sample : db.rawbed_samples )
        {
            sample.second = rawbed_samples_map[ sample.first ];
        }

        monitor.log( "Annotating", " ... Complete" );
        monitor.log( "Component Annotator", "Complete!!" );
    }
};

} // end of namespace component
} // end of namespace ago
