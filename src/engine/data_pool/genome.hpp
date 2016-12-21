#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CCD/utility/language.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <CPT/logger.hpp>
#include <pokemon/iohandler/ihandler/ihandler.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>

namespace ago {
namespace engine {
namespace data_pool {

class GenomeImpl
{
  protected:
    
    bool is_genome_load_;

  public:

    std::map< std::string, std::string > genome_table;

    template< class DB >
    GenomeImpl( DB& db )
    : is_genome_load_( false )
    {
        auto& pipeline_schema ( db.pipeline_schema() );
        for ( auto& child 
                : pipeline_schema.get_child( "input" ).get_child( "genome_fasta" ))
        {
            db.push_path( "genome_fasta", child.second );
        }
    }

    template< class DB >
    void load_genome( std::vector< std::string >& genome_fastas, DB& db )
    {
        auto& monitor = db.monitor();

        monitor.set_monitor( "Loading Genome", genome_fastas.size()+2 );
        monitor.log( "Loading Genome", "Start" );

        std::mutex fa_mutex;
        ParaThreadPool fa_parallel_pool( genome_fastas.size() );

        for( size_t i = 0; i < genome_fastas.size(); ++i )
        {
            fa_parallel_pool.job_post( [ i, &genome_fastas, &db, &fa_mutex, &monitor ] () mutable 
            {
                Fasta_ihandler_impl< IoHandlerIfstream > fasta_reader( genome_fastas );
                Fasta<> chr = fasta_reader.get_next_entry( i );

                {
                    std::lock_guard< std::mutex > fa_lock( fa_mutex );
                    db.genome_table.emplace( chr.getName(), chr.getSeq() );
                    monitor.log( "Loading Genome", ( chr.getName() ).c_str() );
                }
            });
        }

        fa_parallel_pool.flush_pool();
        monitor.log( "Loading Genome", "Complete" );
    }

    template< class DB >
    std::vector< std::string > require_genome( DB& db )
    {
        std::vector< std::string > genome_fastas;

        if( db.exist_path_tag( "genome_fasta" ))
        {
            for( auto& fasta : db.get_path_list( "genome_fasta" ))
            {
                genome_fastas.emplace_back( fasta.string() );
            }
        }
        else
        {
            throw std::runtime_error( "\"genome_fasta\" is not in the db pool" );
        }

        if( !db.is_genome_load_ )
        {
            db.load_genome( genome_fastas, db );
            db.is_genome_load_ = true;
        }
        
        return genome_fastas;
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
