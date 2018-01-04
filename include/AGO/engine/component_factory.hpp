#pragma once
#include <AGO/engine/component.hpp>
#include <memory>
#include <string>
#include <AGO/engine/data_pool.hpp>
#include <AGO/components/tailor_fastq_to_bed.hpp>
#include <AGO/components/tailor_fasta_to_bed.hpp>
#include <AGO/components/fastq_input.hpp>
#include <AGO/components/fasta_input.hpp>
#include <AGO/components/sam_input.hpp>
#include <AGO/components/archive_input.hpp>
#include <AGO/components/tailor_align.hpp>
#include <AGO/components/sam_to_bed.hpp>
#include <AGO/components/annotator.hpp>
#include <AGO/components/filter.hpp>
#include <AGO/components/analyzer.hpp>
#include <AGO/components/vcf_annotation_analysis.hpp>
#include <AGO/components/vcf_min_annotation_cohorts_analysis.hpp>
#include <AGO/components/filter_analyzer.hpp>
#include <AGO/components/meta_analyzer.hpp>
#include <AGO/components/visualization.hpp>
#include <AGO/components/github_tailor_fastq_to_bed.hpp>
#include <AGO/components/github_tailor_fastq_to_sam_out.hpp>

namespace ago {
namespace engine {

#define ID_MAP_TYPE( ID, TYPE ) \
    else if ( identifier == ID ) \
    { \
        auto tmp = ComponentPtr( \
            new TYPE(data_pool_, schema_node) \
        ); \
        tmp->config(schema_node); \
        return tmp; \
    }

class ComponentFactory
{
  public : 

    const DataPool& data_pool_;

    ComponentFactory( const DataPool& data_pool )
    : data_pool_( data_pool )
    {}

    ComponentPtr create_by_identifier ( 
          const std::string& identifier 
        , const bpt::ptree& schema_node
    )
    {

        if( identifier == "" ){ /* TODO a error handle */ }
        ID_MAP_TYPE( "TailorFastqToBed" , component::TailorFastqToBed )
        ID_MAP_TYPE( "TailorFastaToBed" , component::TailorFastaToBed )
        ID_MAP_TYPE( "FastqInput" , component::FastqInput )
        ID_MAP_TYPE( "FastaInput" , component::FastaInput )
        ID_MAP_TYPE( "SamInput" , component::SamInput )
        ID_MAP_TYPE( "ArchiveInput" , component::ArchiveInput )
        ID_MAP_TYPE( "TailorAlign" , component::TailorAlign )
        ID_MAP_TYPE( "SamToBed" , component::SamToBed )
        ID_MAP_TYPE( "Annotator" , component::Annotator )
        ID_MAP_TYPE( "Filter" , component::Filter )
        ID_MAP_TYPE( "Analyzer", component::Analyzer )
        ID_MAP_TYPE( "VcfAnnotationAnalysis", component::VcfAnnotationAnalysis )
        ID_MAP_TYPE( "VcfMinAnnotationCohortsAnalysis", component::VcfMinAnnotationCohortsAnalysis )
        ID_MAP_TYPE( "FilterAnalyzer", component::FilterAnalyzer )
        ID_MAP_TYPE( "MetaAnalyzer", component::MetaAnalyzer )
        ID_MAP_TYPE( "Visualization", component::Visualization )
        ID_MAP_TYPE( "GithubTailorFastqToBed", component::GithubTailorFastqToBed )
        ID_MAP_TYPE( "GithubTailorFastqToSamOut", component::GithubTailorFastqToSamOut )
        /* TODO ADD COMPONENT HERE */
        else { /* TODO another error handle */ }

        return nullptr;
    }

    template<class T>
    auto operator () ( const T& component_schema )
    {
        auto name = component_schema
            .second
            .template get<std::string> ("name"); 
        return create_by_identifier ( 
              name
            , component_schema.second
        );
    }
};

#undef ID_MAP_TYPE

} // engine
} // ago
