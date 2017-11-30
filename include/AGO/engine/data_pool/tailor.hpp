#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class Tailor
{
  public:

    bool is_tailor_index_build;

    template< class DB >
    Tailor( DB& db )
    : is_tailor_index_build( false )
    {
    }

    template< class DB >
    std::string require_tailor_index( const std::vector< std::string >& genome_fastas, DB& db )
    {
        std::string tailor_index;

        if( db.exist_path_tag( "tailor_index" ))
        {
            tailor_index = db.get_path( "tailor_index" ).string();
        }
        else
        {
            throw std::runtime_error( "\"tailor_index\" is not in the db pool" );
        }

        if( boost::filesystem::exists(
                boost::filesystem::path( tailor_index + ".t_table.bwt" )
        ))
        {
            db.is_tailor_index_build = true;
        }

        return tailor_index;
    }

    template< class DB >
    std::string require_tailor_index( DB& db )
    {
        std::string tailor_index;

        if( db.exist_path_tag( "tailor_index" ))
        {
            tailor_index = db.get_path( "tailor_index" ).string();
        }
        else
        {
            throw std::runtime_error( "\"tailor_index\" is not in the db pool" );
        }

        if( boost::filesystem::exists(
                boost::filesystem::path( tailor_index + "t_table.bwt" )
        ))
        {
            db.is_tailor_index_build = true;
        }

        return tailor_index;
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
