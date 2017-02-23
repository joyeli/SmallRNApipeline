#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class FastaInput : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t reads_min_length_;
    size_t reads_max_length_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            db.push_path( "sample_files", child.second );
        }

        reads_min_length_ = p.get_optional< size_t >( "min_length" ).value_or( 15 );
        reads_max_length_ = p.get_optional< size_t >( "max_length" ).value_or( 30 );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > fasta_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fasta_ihandler_impl< IoHandlerIfstream > fasta_reader( fasta_paths );

        monitor.set_monitor( "Component FastaInput", fasta_paths.size()+2 );
        monitor.log( "Component FastaInput", "Start" );

        std::mutex fq_mutex;
        ParaThreadPool fq_parallel_pool( fasta_paths.size() );

        for( auto& fasta_path : fasta_paths )
        {
            std::string sample_name( get_sample_name( fasta_path ));
            db.fastq_samples.emplace_back( sample_name, std::vector< Fastq<> >{} );
        }

        for( size_t smp = 0; smp < fasta_paths.size(); ++smp )
        {
            fq_parallel_pool.job_post( [ smp, &db, &fasta_paths, &fasta_reader, &fq_mutex, &monitor, this ] ()
            {
                std::vector< Fastq<> > fastqs;
                std::string qc = "";
                Fastq<> fastq;
                Fasta<> fasta;

                while( true )
                {
                    fasta = fasta_reader.get_next_entry( smp );

                    if( fasta.eof_flag )
                    {
                        break;
                    }

                    if( fasta.getSeq().size() < reads_min_length_ ||
                        fasta.getSeq().size() > reads_max_length_ ||
                        n_check( fasta )
                    )
                    {
                        continue;
                    }

                    for( auto& c : std::get<1>( fasta.data ))
                    {
                        qc += "I";
                    }

                    switch( std::get<0>( fasta.data ).at( 0 ))
                    {
                        case '@':
                            std::get<0>( fastq.data ) = std::get<0>( fasta.data );
                            break;

                        default :
                            std::get<0>( fastq.data ) = "@" + std::get<0>( fasta.data );
                            break;
                    }

                    std::get<1>( fastq.data ) = std::get<1>( fasta.data );
                    std::get<2>( fastq.data ) = "+" + std::get<0>( fasta.data );
                    std::get<3>( fastq.data ) = qc;
                    qc = "";

                    fastqs.push_back( fastq );
                }

                {
                    std::lock_guard< std::mutex > fq_lock( fq_mutex );
                    db.fastq_samples[ smp ].second = std::move( fastqs );
                    monitor.log( "Component FastaInput", ( db.fastq_samples[ smp ].first ).c_str() );
                }
            });
        }

        fq_parallel_pool.flush_pool();
        monitor.log( "Component FastaInput", "Complete" );
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

    bool n_check( Fasta<>& fasta )
    {
        for( auto& base : fasta.getSeq() )
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
