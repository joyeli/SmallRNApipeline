#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <pokemon/converter/sam2rawbed.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>

namespace ago {
namespace component {

class TailorFastaToBed : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    size_t reads_min_length_;
    size_t reads_max_length_;

    size_t align_min_length_;
    size_t align_max_length_;
    size_t align_limit_algn_;

    size_t task_number_;
    size_t thread_num_;

    Aligner< Aligner_trait< Aligner_types::BWT_Aligner, ParallelTypes::M_T, void >> tailor_;

    bool output_sam_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
        {
            db.push_path( "sample_files", child.second );
        }

        std::string tailor_index = p.get_optional< std::string >( "tailor_index" ).value_or( "" );

        if( tailor_index != "" && !db.exist_path_tag( "tailor_index" ))
        {
            db.push_path( "tailor_index", tailor_index );
        }

        reads_min_length_ = p.get_optional< size_t >( "reads_min_length" ).value_or( 15 );
        reads_max_length_ = p.get_optional< size_t >( "reads_max_length" ).value_or( 30 );

        align_min_length_ = p.get_optional< size_t >( "align_min_length" ).value_or( 18 );
        align_max_length_ = p.get_optional< size_t >( "align_max_length" ).value_or( 30 );
        align_limit_algn_ = p.get_optional< size_t >( "align_limit_algn" ).value_or( 100 );

        task_number_ = p.get_optional< size_t >( "task_number" ).value_or( 1000000 );
        thread_num_ = p.get_optional< size_t >( "thread_num" ).value_or( 32 );

        output_sam_ = p.get_optional< bool >( "output_sam" ).value_or( false );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        if( !db.exist_path_tag( "tailor_index" ))
        {
            auto path = db.output_dir();
            path += "tailor";
            db.push_path( "tailor_index", path );
        }

        std::vector< std::string > genome_fastas( db.require_genome( db ));
        std::string tailor_index( db.require_tailor_index( genome_fastas, db ));

        if( db.is_tailor_index_build )
        {
            monitor.set_monitor( "Loading Index", 2 );
            monitor.log( "Loading Index", "Start" );

            tailor_.load_table( tailor_index );
            monitor.log( "Loading Index", "Complete" );
        }
        else
        {
            monitor.set_monitor( "Building Index", 3 );
            monitor.log( "Building Index", "Building ... " );

            tailor_.build( genome_fastas, tailor_index );
            monitor.log( "Building Index", "Loading ... " );

            tailor_.load_table( tailor_index );
            db.is_tailor_index_build = true;

            monitor.log( "Building Index", "Complete" );
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        std::vector< std::string > fasta_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fasta_ihandler_impl< IoHandlerIfstream > fasta_reader( fasta_paths );

        monitor.set_monitor( "Component TailorFastaToBed", fasta_paths.size() +2 );
        monitor.log( "Component TailorFastaToBed", "Start" );

        std::mutex smp_mutex;
        ParaThreadPool smp_parallel_pool( fasta_paths.size() );


        for( size_t id = 0; id < fasta_paths.size(); ++id )
        {
            db.statistic_samples.emplace_back( "", std::vector< double >{ 0.0, 0.0, 0.0 });

            smp_parallel_pool.job_post( [ id, &db, &fasta_paths, &fasta_reader, &smp_mutex, &monitor, this ] ()
            {
                std::string sample_name( get_sample_name( fasta_paths[ id ] ));

                std::mutex ali_mutex;
                ParaThreadPool ali_parallel_pool( thread_num_ );
                
                std::ofstream output_sam;

                if( output_sam_ )
                {
                    output_sam.open( db.output_dir().string() + sample_name + ".sam" );
                }

                std::vector< Sam<> > sams{};

                size_t thread_count = 0;
                bool break_flag = false;

                std::vector< Fastq<> > fastqs;
                std::string qc = "";
                Fastq<> fastq;
                Fasta<> fasta;

                std::map< std::string, size_t > align_count;
                std::map< std::string, size_t > fastq_count;
                std::map< std::string, size_t >::iterator fq_it;

                while( true )
                {
                    for( size_t job = 0; job < task_number_; ++job )
                    {
                        fasta = fasta_reader.get_next_entry( id );

                        if( fasta.eof_flag )
                        {
                            break_flag = true;
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

                        std::get<0>( fastq.data ) = std::get<0>( fasta.data );
                        std::get<1>( fastq.data ) = std::get<1>( fasta.data );
                        std::get<2>( fastq.data ) = "+" + std::get<0>( fasta.data );
                        std::get<3>( fastq.data ) = qc;
                        qc = "";

                        fq_it = fastq_count.find( std::get<0>( fastq.data ));

                        if( fq_it != fastq_count.end() )
                        {
                            fq_it->second++;
                        }
                        else
                        {
                            fastq_count.emplace( std::get<0>( fastq.data ), 1 );
                        }

                        fastqs.push_back( fastq );
                    }

                    ali_parallel_pool.job_post( [ fastqs, &ali_mutex, &sams, &output_sam, this ] () 
                    {
                        std::vector< Sam<> > sams_para{};
                        std::map< int, std::vector< Fastq<> >> input;
                        std::vector< Sam<> >* sams_tmp;

                        for( auto& fastq: fastqs )
                        {
                            input.emplace( 0, std::vector< Fastq<> >{ fastq });
                            sams_tmp = tailor_.search( &input, 1, align_min_length_, align_limit_algn_, align_max_length_ );

                            if( !sams_tmp->empty() )
                            {
                                std::move( sams_tmp->begin(), sams_tmp->end(), std::back_inserter( sams_para ));
                            }

                            input.clear();
                            sams_tmp->clear();
                        }

                        {
                            std::lock_guard< std::mutex > ali_lock( ali_mutex );

                            if( output_sam_ )
                            {
                                for( auto& sam : sams_para )
                                {
                                    output_sam << sam;
                                }
                            }

                            std::move( sams_para.begin(), sams_para.end(), std::back_inserter( sams ));
                        }
                    });

                    thread_count++;
                    fastqs.clear();

                    if( thread_count >= thread_num_ )
                    {
                        thread_count = 0;
                        ali_parallel_pool.flush_pool();
                    }

                    if( break_flag )
                    {
                        ali_parallel_pool.flush_pool();

                        if( output_sam_ )
                        {
                            output_sam.close();
                        }

                        break;
                    }
                }

                for( auto sam = sams.begin(); sam != sams.end(); ++sam )
                {
                    std::get<3>( sam->data ) = std::get<3>( sam->data ) - 1;
                    align_count.emplace( "@" + std::get<0>( sam->data ), 0 );
                }

                {
                    std::lock_guard< std::mutex > smp_lock( smp_mutex );

                    Sam2RawBed< std::vector< Sam<> >* > sam2bed;
                    auto raw_beds( sam2bed.Convert( &sams ));

                    std::vector< AnnotationRawBed<> > annotation_rawbeds;

                    for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
                    {
                        annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
                    }

                    sam2bed.rawbed_map_->clear();
                    sam2bed.rawbed_map2_->clear();

                    db.bed_samples.emplace_back( sample_name, annotation_rawbeds );

                    db.statistic_samples[ id ].first = sample_name;

                    for( auto& fqc : fastq_count )
                    {
                        db.statistic_samples[ id ].second[0] += fqc.second;

                        if( align_count.find( fqc.first ) != align_count.end() )
                        {
                            db.statistic_samples[ id ].second[1] += fqc.second;
                        }
                    }

                    db.statistic_samples[ id ].second[2] =
                        db.statistic_samples[ id ].second[1] * 100 / db.statistic_samples[ id ].second[0];

                    monitor.log( "Component TailorFastaToBed", ( sample_name ).c_str() );
                }
            });
        }

        smp_parallel_pool.flush_pool();
        make_statistic( db );

        monitor.log( "Component TailorFastaToBed", "Complete" );
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

    void make_statistic( auto& db )
    {
        std::ofstream output( db.output_dir().string() + "mappability.tsv" );

        output << "Sample";

        for( auto& sts : db.statistic_samples )
        {
            output << "\t" << sts.first;
        }

        output << "\n";

        for( size_t i = 0; i < 3; ++i )
        {
            switch( i )
            {
                case 0 : output << "RawRead:"; break;
                case 1 : output << "Mappable:"; break;
                case 2 : output << "Mappable%:"; break;
                default: std::runtime_error( "out of row in mappability.tsv" );
            }

            for( auto& sts : db.statistic_samples )
            {
                output << "\t" << sts.second[i];
            }
            
            output << "\n";
        }

        output.close();
    }
};

} // end of namespace component
} // end of namespace ago
