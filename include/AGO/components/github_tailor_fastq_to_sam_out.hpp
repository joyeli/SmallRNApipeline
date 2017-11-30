#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <Tailor/tailer.hpp>
#include <pokemon/converter/sam2rawbed.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>
#include <set>

#include <sys/resource.h>

namespace ago {
namespace component {

class GithubTailorFastqToSamOut : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    int reads_min_length_;
    int reads_max_length_;

    int align_min_length_;
    int align_limit_algn_;

    int task_number_;
    int thread_num_;

    bool align_allow_mismatch_;

    std::stringstream sam_header_ss_;
    std::string tailor_genome_fasta_;
    ABWT_table abwtt_;

    struct rusage r_usage;
    std::mutex mem_mutex;

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

        reads_min_length_ = p.get_optional< int >( "reads_min_length" ).value_or( 15 );
        reads_max_length_ = p.get_optional< int >( "reads_max_length" ).value_or( 30 );
        align_min_length_ = p.get_optional< int >( "align_min_length" ).value_or( 12 );
        align_limit_algn_ = p.get_optional< int >( "align_limit_algn" ).value_or( 10 );

        task_number_ = p.get_optional< int >( "task_number" ).value_or( 50000 );
        thread_num_  = p.get_optional< int >( "thread_num" ).value_or( 16 );

        align_allow_mismatch_ = p.get_optional< bool >( "align_allow_mismatch" ).value_or( false );
        tailor_genome_fasta_  = p.get_optional< std::string >( "tailor_genome_fasta" ).value_or( "" );
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

        // print_mem_usage( "Initialize Start" );

        std::vector< std::string > genome_fastas( db.require_genome( db ));
        std::string tailor_index( db.require_tailor_index( db ));

        // print_mem_usage( "Genome Fasta Loading" );

        if( db.is_tailor_index_build )
        {
            monitor.set_monitor( "Loading Index", 2 );
            monitor.log( "Loading Index", "Start" );

            abwtt_ = tailor::loadBWT2( tailor_index, &sam_header_ss_ );
            monitor.log( "Loading Index", "Complete" );
        }
        else
        {

            monitor.set_monitor( "Building Index", 3 );
            monitor.log( "Building Index", "Building ..." );

            if( tailor_genome_fasta_ == "" )
                throw std::runtime_error( "\"tailor_genome_fasta\" is required as single fasta file" );

            tailor::buildBWT2( tailor_genome_fasta_, tailor_index );
            monitor.log( "Building Index", "Loading ..." );

            abwtt_ = tailor::loadBWT2( tailor_index, &sam_header_ss_ );
            db.is_tailor_index_build = true;

            monitor.log( "Building Index", "Complete" );
        }

        // print_mem_usage( "Tailor Index Loading" );
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        size_t load_counts = 0;
        // print_mem_usage( "Component Start" );

        std::vector< std::string > fastq_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        Fastq_ihandler_impl< IoHandlerIfstream > fastq_reader( fastq_paths );

        monitor.set_monitor( "Component GithubTailorFastqToSamOut", fastq_paths.size() +2 );
        monitor.log( "Component GithubTailorFastqToSamOut", "Start" );

        sam_header_ss_.clear();
        sam_header_ss_ << "@HD" << '\t' << "VN:1.0" << '\t' << "SO:unsorted\n";

        std::ofstream output_sam;

        std::string line;
        std::string sample_name;

        std::stringstream fastq_ss;

        std::vector< std::string > id_split;
        std::set< std::string > align_count;

        std::map< std::string, size_t > fastq_count;
        std::map< std::string, size_t >::iterator fq_it;

        bool break_flag = false;
        bool new_line_check = false;

        std::mutex sam_mutex;
        ParaThreadPool parallel_pool( thread_num_ );

        for( size_t smp = 0; smp < fastq_paths.size(); ++smp )
        {
            load_counts = 0;

            sample_name = get_sample_name( fastq_paths[ smp ] );
            db.statistic_samples.emplace_back( sample_name, std::vector< double >{ 0.0, 0.0, 0.0 });

            // print_mem_usage( "Sample " + sample_name + " Start" );

            monitor.set_monitor(( "	" + sample_name ).c_str(), 2 );
            monitor.log(( "	" + sample_name ).c_str(), "Aligning ..." );

            output_sam.open( db.output_dir().string() + sample_name + ".sam" );
            output_sam << sam_header_ss_.str();

            break_flag = false;

            while( true )
            {
                load_counts++;
                new_line_check = false;

                for( size_t job = 0; job < task_number_; ++job )
                {
                    Fastq<> fastq( fastq_reader.get_next_entry( smp ));

                    if( fastq.eof_flag )
                    {
                        break_flag = true;
                        break;
                    }

                    if( fastq.getSeq().size() < reads_min_length_ ||
                        fastq.getSeq().size() > reads_max_length_ ||
                        n_check( fastq )
                      )
                    {
                        --job;
                        continue;
                    }

                    fq_it = fastq_count.find( fastq.getName() );

                    if( fq_it != fastq_count.end() )
                    {
                        fq_it->second++;
                    }
                    else
                    {
                        line = fastq.getName();
                        boost::iter_split( id_split, line, boost::algorithm::first_finder( " " ));
                        fastq_count.emplace( id_split[0], 1 );
                    }

                    if( new_line_check )
                    {
                        fastq_ss << "\n";
                        fastq_ss.clear();
                    }
                    else
                    {
                        new_line_check = true;
                    }

                    fastq_ss
                        << "@" << fastq.getName() << "\n"
                        << fastq.getSeq() << "\n"
                        << fastq.getName2() << "\n"
                        << fastq.getQuality();
                }

                if( !new_line_check && break_flag )
                {
                    break;
                }

                // print_mem_usage( "Sample " + sample_name + " Fastq " + std::to_string( load_counts * task_number_ ) + " Loading-" + std::to_string( load_counts ));

                parallel_pool.job_post([ fastq_str = fastq_ss.str(), &sam_mutex, load_counts, &output_sam, &align_count, this ] () mutable
                {
                    aligning( fastq_str, sam_mutex, load_counts, output_sam, align_count );
                });

                fastq_ss.str("");
                fastq_ss.clear();

                // print_mem_usage( "Sample " + sample_name + " Fastq " + std::to_string( load_counts * task_number_ ) + " Releasing-" + std::to_string( load_counts ));

                if( break_flag )
                {
                    break;
                }
            }

            parallel_pool.flush_pool();
            output_sam.close();

            // print_mem_usage( "Sample " + sample_name + " Aligning" );

            for( auto& fqc : fastq_count )
            {
                db.statistic_samples[ smp ].second[0] += fqc.second;

                if( align_count.find( fqc.first ) != align_count.end() )
                {
                    db.statistic_samples[ smp ].second[1] += fqc.second;
                }
            }

            fastq_count.clear();
            align_count.clear();

            db.statistic_samples[ smp ].second[2] =
                db.statistic_samples[ smp ].second[1] * 100 / db.statistic_samples[ smp ].second[0];

            // print_mem_usage( "Sample " + sample_name + " Statistic Releasing" );

            monitor.log(( "	" + sample_name ).c_str(), "Done" );
            monitor.log( "Component GithubTailorFastqToSamOut", ( sample_name ).c_str() );

            // print_mem_usage( "Sample " + sample_name + " Complete" );
        }

        make_statistic( db );
        monitor.log( "Component GithubTailorFastqToSamOut", "Complete" );

        // print_mem_usage( "Component Complete" );
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
                output << std::fixed << std::setprecision( 0 ) << "\t" << sts.second[i];
            }
            
            output << "\n";
        }

        output.close();
    }

    void aligning( auto& fastq_str, auto& sam_mutex, auto& load_counts, auto& output_sam, auto& align_count )
    {
        // print_mem_usage( "Tailor Start-" + std::to_string( load_counts ));

        std::stringstream sam_ss;
        std::stringstream fastq_ss;
        boost::thread* threads[1];

        fastq_ss << fastq_str;
        fastq_str = "";

        threads[0] = new boost::thread
        {
            tailor::ABWT_threads< ABWT_table >
            {
                abwtt_, &fastq_ss, &sam_ss, align_min_length_, align_allow_mismatch_ 
            }
        };

        if( threads[0]->joinable() )
        {
            threads[0]->join();
        }

        // print_mem_usage( "Tailor Aligning-" + std::to_string( load_counts ));

        delete threads[0];
        fastq_ss.str("");
        fastq_ss.clear();

        // print_mem_usage( "Tailor Fastq Releasing-" + std::to_string( load_counts ));

        std::string samline;
        SamDefaultTuple sam_tpl;
        Sam< SamDefaultTuple > sam;

        std::vector< std::string > split;
        std::vector< std::string > sam_nh;
        std::vector< std::pair< std::string, std::string >> sams;

        while( true )
        {
            if( !std::getline( sam_ss, samline ))
            {
                break;
            }

            boost::iter_split( split, samline, boost::algorithm::first_finder( "\t" ));
            boost::iter_split( sam_nh, split[11], boost::algorithm::first_finder( ":" ));

            if( std::stoi( sam_nh[2] ) > align_limit_algn_ )
            {
                continue;
            }

            sams.emplace_back( std::make_pair( split[0], samline ));

            split.clear();
            sam_nh.clear();
        }

        // print_mem_usage( "Tailor Sam Converting-" + std::to_string( load_counts ));

        sam_ss.str("");
        sam_ss.clear();

        for( auto& sam : sams )
        {
            std::lock_guard< std::mutex > sam_lock( sam_mutex );
            align_count.emplace( sam.first );
            output_sam << sam.second << "\n";
        }

        sams.clear();

        // print_mem_usage( "Tailor Sam Releasing-" + std::to_string( load_counts ));
    }

    void print_mem_usage( std::string msg )
    {
        std::lock_guard< std::mutex > mem_lock( mem_mutex );

        double vm_usage     = 0.0;
        double resident_set = 0.0;
        // 'file' stat seems to give the most reliable results
        std::ifstream stat_stream( "/proc/self/stat", std::ios_base::in );

        // dummy vars for leading entries in stat that we don't care about
        std::string pid, comm, state, ppid, pgrp, session, tty_nr;
        std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
        std::string utime, stime, cutime, cstime, priority, nice;
        std::string O, itrealvalue, starttime;

        // the two fields we want
        size_t vsize;
        size_t rss;

        stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
            >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
            >> utime >> stime >> cutime >> cstime >> priority >> nice
            >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

        stat_stream.close();

        size_t page_size_kb = sysconf( _SC_PAGE_SIZE ); // in case x86-64 is configured to use 2MB pages

        vm_usage     = vsize / 1073741824.0;
        resident_set = rss * page_size_kb / 1073741824.0;

        getrusage( RUSAGE_SELF, &r_usage );
        std::cout << msg << " ... ... ... ... VM: " << vm_usage << " GB / RSS: " << resident_set << " GB / HIGHTST: " << double( r_usage.ru_maxrss )/ 1048576.0 << " GB" << std::endl;
    }
};

} // end of namespace component
} // end of namespace ago
