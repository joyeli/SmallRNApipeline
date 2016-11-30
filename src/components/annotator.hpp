#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>

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
            if( !db.exist_path_tag( "annotation_files" ))
            {
                db.push_path( "annotation_files", child.second );
            }
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

        Annotations annotator( annotation_files_ );

        for( auto& sample : db.rawbed_samples )
        {
            for( auto& anno_rawbed : sample.second )
            {
                annotator.AnnotateAll( anno_rawbed );
            }
        }

        Annotations::clear_database();
    }
};

} // end of namespace component
} // end of namespace ago
