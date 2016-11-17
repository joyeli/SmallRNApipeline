#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <CCD/utility/language.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <CPT/logger.hpp>
#include <utility>
#include <set>
#include <pokemon/iohandler/ihandler/ihandler.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class FastqImpl
{
  protected:

    std::set< std::string > loaded_fastq_;

  public:

    std::vector<
        std::pair<
            std::string, // sample name
            std::vector< // fastq vector
                Fastq<>
            >
        >
    > fastq_samples; 

    size_t fastq_min_length;
    size_t fastq_max_length;

    template< class DB >
    FastqImpl( DB& db )
    {
        auto& pipeline_schema ( db.pipeline_schema() );
        for ( auto& child 
                : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            if( !db.exist_path_tag( "sample_files" ))
            {
                db.push_path( "sample_files", child.second );
            }
        }
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

    template< class DB >
    void load_fastq( const std::string& path, DB& db )
    {
        cpt::verbose0 << "load " << path << std::endl;

        std::vector< Fastq<> > fastqs;
        std::vector< std::string > f_path{ path };

        Fastq_ihandler_impl< IoHandlerIfstream > fastq_reader( f_path );

        while( true )
        {
            Fastq<> fastq = fastq_reader.get_next_entry( 0 );

            if( fastq.eof_flag )
            {
                break;
            }

			if( fastq.getSeq().size() < fastq_min_length ||
                fastq.getSeq().size() > fastq_max_length ||
                n_check( fastq )
            )
				continue;

            fastqs.push_back( fastq );
        }

        std::vector< std::string > path_file;
        boost::iter_split( path_file, path, boost::algorithm::first_finder( "/" ));

        std::vector< std::string > fastq_file;
        boost::iter_split( fastq_file, path_file[ path_file.size()-1 ], boost::algorithm::first_finder( "." ));

        db.fastq_samples.emplace_back( std::make_pair( fastq_file[0], fastqs )); 
    }

    template< class DB >
    void require_fastq( const std::string& path, DB& db )
    {
        if( db.loaded_fastq_.find( path ) == db.loaded_fastq_.end() )
        {
            cpt::verbose0 << "load : " << path << std::endl;
            db.load_fastq( path, db );
            db.loaded_fastq_.emplace( path );
        }
    }
};


} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
