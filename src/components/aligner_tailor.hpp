#pragma once
#include <AGO/engine/components/named_component.hpp>

namespace ago {
namespace component {

class AlignerTailor : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    typedef Aligner_trait< 
        Aligner_types::BWT_Aligner, 
        ParallelTypes::M_T,
        void
    >  trait;

    Aligner< trait > tailor_;

    size_t align_thread_num_;
    size_t align_min_length_;
    size_t align_limit_num_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );

        align_thread_num_ = p.get_optional< size_t >( "align_thread_num" ).value_or( 10000 );
        align_min_length_ = p.get_optional< size_t >( "align_min_length" ).value_or( 18 );
        align_limit_num_ = p.get_optional< size_t >( "align_min_length" ).value_or( 100 );

        std::string tailor_index = p.get_optional< std::string >( "tailor_index" ).value_or( "" );

        if( tailor_index != "" && !db.exist_path_tag( "tailor_index" ))
        {
            db.push_path( "tailor_index", tailor_index );
        }
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );

        if( !db.exist_path_tag( "tailor_index" ))
        {
            auto path = db.output_dir();
            path += "/tailor";
            db.push_path( "tailor_index", path );
        }

        std::string tailor_index( db.require_tailor_index( db ));
        std::vector< std::string > genome_fasta( db.require_genome( db ));

        if( !db.is_tailor_index_build )
        {
            cpt::verbose0 << "build index : " << tailor_index << std::endl;

            tailor_.build( genome_fasta, tailor_index );

            db.is_tailor_index_build = true;
        }
        else
        {
            cpt::verbose0 << "load index : " << tailor_index + ".t_table.bwt" << std::endl;

            tailor_.load_table( tailor_index );
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );

        for( auto& fastq_sample : db.fastq_samples )
        {
            std::map< int, std::vector< Fastq<> >> fastq_maps{ std::make_pair( 0, std::move( fastq_sample.second ))};
            std::vector< Sam<> >* sams( tailor_.search( &fastq_maps, align_thread_num_, align_min_length_, align_min_length_ ));

            db.sam_samples.emplace_back( fastq_sample.first, sams );
        }
    }
};

} // end of namespace component
} // end of namespace ago
