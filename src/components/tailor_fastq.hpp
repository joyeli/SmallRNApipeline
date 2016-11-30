#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <pokemon/converter/sam2rawbed.hpp>

namespace ago {
namespace component {

class TailorFastq : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t fastq_min_length_;
    size_t fastq_max_length_;

    size_t align_min_length_;
    size_t align_max_length_;

    size_t align_job_number_;
    size_t align_thread_num_;
    size_t align_limit_algn_;

    Aligner< Aligner_trait< Aligner_types::BWT_Aligner, ParallelTypes::M_T, void >> tailor_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            if( !db.exist_path_tag( "sample_files" ))
            {
                db.push_path( "sample_files", child.second );
            }
        }

        std::string tailor_index = p.get_optional< std::string >( "tailor_index" ).value_or( "" );

        if( tailor_index != "" && !db.exist_path_tag( "tailor_index" ))
        {
            db.push_path( "tailor_index", tailor_index );
        }

        fastq_min_length_ = p.get_optional< size_t >( "fastq_min_length" ).value_or( 15 );
        fastq_max_length_ = p.get_optional< size_t >( "fastq_max_length" ).value_or( 30 );

        align_min_length_ = p.get_optional< size_t >( "align_min_length" ).value_or( 18 );
        align_max_length_ = p.get_optional< size_t >( "align_max_length" ).value_or( 30 );

        align_job_number_ = p.get_optional< size_t >( "align_job_number" ).value_or( 20000 );
        align_thread_num_ = p.get_optional< size_t >( "align_thread_num" ).value_or( 200 );
        align_limit_algn_ = p.get_optional< size_t >( "align_limit_algn" ).value_or( 100 );
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

        std::vector< std::string > genome_fastas( db.require_genome( db ));
        std::string tailor_index( db.require_tailor_index( genome_fastas, db ));

        if( db.is_tailor_index_build )
        {
            cpt::verbose0 << "load index : " << tailor_index + ".t_table.bwt" << std::endl;
            tailor_.load_table( tailor_index );
        }
        else
        {
            cpt::verbose0 << "build index : " << tailor_index << std::endl;
            tailor_.build( genome_fastas, tailor_index );
            db.is_tailor_index_build = true;
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );

        std::vector< std::string > fastq_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fastq_ihandler_impl< IoHandlerIfstream > fastq_reader( fastq_paths );

        for( size_t id = 0; id < fastq_paths.size(); ++id )
        {
            cpt::verbose0 << "load : " << fastq_paths[ id ] << std::endl;
            std::string sample_name( get_sample_name( fastq_paths[ id ] ));

            std::vector< Sam<> >* sams;

            bool break_flag = false;

            while( true )
            {
                std::vector< Fastq<> > fastqs;

                for( size_t job = 0; job < align_job_number_; ++job )
                {
                    Fastq<> fastq = fastq_reader.get_next_entry( id );

                    if( fastq.eof_flag )
                    {
                        break_flag = true;
                        break;
                    }

                    if( fastq.getSeq().size() < fastq_min_length_ ||
                        fastq.getSeq().size() > fastq_max_length_ ||
                        n_check( fastq )
                      )
                        continue;

                    fastqs.push_back( fastq );
                }

                std::map< int, std::vector< Fastq<> >> input{ std::make_pair( 0, fastqs )};

                std::vector< Sam<> >* sams_tmp( tailor_.search( &input, align_thread_num_, align_min_length_, align_limit_algn_, align_max_length_ ));
                sams->insert(  sams->end(), sams_tmp->begin(), sams_tmp->end() );

                fastqs.clear();

                if( break_flag )
                {
                    break;
                }
            }

            auto raw_beds( Sam2RawBed< std::vector< Sam<> >* >().Convert( sams ));

            std::vector< AnnotationRawBed<> > annotation_rawbeds;

            for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
            {
                annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
            }

            db.rawbed_samples.emplace_back( sample_name, annotation_rawbeds );
        }
    }

    std::vector< std::string > get_path_list_string( const std::vector< boost::filesystem::path >& paths )
    {
        std::vector< std::string > res;

        for( auto& path : paths )
        {
            res.emplace_back( path.string() );
        }

        return res;
    }

    std::string get_sample_name( const std::string& path )
    {
        std::vector< std::string > path_file;
        boost::iter_split( path_file, path, boost::algorithm::first_finder( "/" ));

        std::vector< std::string > sample;
        boost::iter_split( sample, path_file[ path_file.size()-1 ], boost::algorithm::first_finder( "." ));

        return sample[0];
    }

    bool n_check( const Fastq<>& fastq )
    {
        for( auto& base : fastq.getSeq() )
        {
            switch( base )
            {
                case 'A': break;
                case 'T': break;
                case 'C': break;
                case 'G': break;
                default : return true;
            }
        }

        return false;
    }
};

} // end of namespace component
} // end of namespace ago
