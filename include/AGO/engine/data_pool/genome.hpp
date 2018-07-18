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
    void load_genome_fasta( std::vector< std::string >& genome_fastas, std::string& genome_arc, DB& db )
    {
        auto& monitor = db.monitor();

        std::ofstream archive_output( genome_arc );
        boost::archive::binary_oarchive archive_out( archive_output );

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

        archive_out & db.genome_table;
        archive_output.close();

        monitor.log( "Loading Genome", "Complete" );
    }

    template< class DB >
    void load_genome_arc( std::string& genome_arc, DB& db )
    {
        auto& monitor = db.monitor();

        monitor.set_monitor( "Loading Genome", 2 );
        monitor.log( "Loading Genome", "Start" );

        std::ifstream archive( genome_arc );
        boost::archive::binary_iarchive archive_in( archive );

        archive_in & db.genome_table;
        archive.close();

        monitor.log( "Loading Genome", "Complete" );
    }

    template< class DB >
    void require_genome( DB& db )
    {
        bool is_archive = false;
        std::string genome_arc = "";

        std::vector< std::string > split;
        std::vector< std::string > genome_fastas;

        if( db.exist_path_tag( "genome_fasta" ))
        {
            for( auto& fasta : db.get_path_list( "genome_fasta" ))
            {
                genome_fastas.emplace_back( fasta.string() );
            }

            boost::iter_split( split, genome_fastas[0], boost::algorithm::first_finder( "/" ));

            for( std::size_t i = 0; i < split.size() -1; ++i )
                genome_arc += split[i] + "/";

            genome_arc += "genomes.arc";

            if( boost::filesystem::exists( genome_arc ))
                is_archive = true;
        }
        else
        {
            throw std::runtime_error( "\"genome_fasta\" or \"genome_arc\" are not in the db pool" );
        }

        if( !db.is_genome_load_ )
        {
            if( is_archive )
                db.load_genome_arc( genome_arc, db );
            else
                db.load_genome_fasta( genome_fastas, genome_arc, db );

            db.is_genome_load_ = true;
        }
    }
};

} // end of namespace data_pool
} // end of namespace engine
} // end of namespace cpt
