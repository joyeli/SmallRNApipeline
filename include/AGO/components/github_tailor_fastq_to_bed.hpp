#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <pokemon/aligner/aligner.hpp>
#include <Tailor/tailer.hpp>
// #include <pokemon/converter/sam2rawbed.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <mutex>
#include <set>

#include <AGO/format/md_sam.hpp>
#include <AGO/format/md_rawbed.hpp>
#include <sys/resource.h>

namespace ago {
namespace component {

class GithubTailorFastqToBed : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

    int align_min_length_;
    int align_min_multi_;
    std::size_t max_tail_len;

    int task_number_;
    int thread_num_;

    bool output_sam_;
    bool align_allow_mismatch_;
    bool align_allow_t2c_;

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

        align_min_length_ = p.get_optional< int >( "align_min_length" ).value_or( 12 );
        align_min_multi_  = p.get_optional< int >( "align_min_multi"  ).value_or( 10 );

        max_tail_len = p.get_optional< std::size_t >( "max_tail_len" ).value_or( 5 );
        task_number_ = p.get_optional< int >( "task_number" ).value_or( 50000 );
        thread_num_  = p.get_optional< int >( "thread_num" ).value_or( 16 );

        output_sam_ = p.get_optional< bool >( "output_sam" ).value_or( true );

        align_allow_mismatch_ = p.get_optional< bool >( "align_allow_mismatch" ).value_or( false );
        align_allow_t2c_      = p.get_optional< bool >( "align_allow_t2c" ).value_or( false );
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

        db.require_genome( db );
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

        monitor.set_monitor( "Component GithubTailorFastqToBed", fastq_paths.size() +2 );
        monitor.log( "Component GithubTailorFastqToBed", "Start" );

        sam_header_ss_.clear();
        sam_header_ss_ << "@HD" << '\t' << "VN:1.0" << '\t' << "SO:unsorted\n";

        std::ofstream output_sam;
        std::stringstream fastq_ss;

        std::string sample_name;
        std::set< std::string > align_count;
        // std::vector< Sam<> > sams{};


        std::map< std::string, size_t > rawread_count;
        std::map< std::string, size_t > fastq_n_count;
        std::map< std::string, size_t > fastq_count;

        bool break_flag = false;
        bool new_line_check = false;
        bool inner_break_flag = false;

        std::mutex sam_mutex;
        ParaThreadPool parallel_pool( thread_num_ );

        std::vector< std::pair< std::string, std::vector< double >>> statistic_samples;
        std::map< ago::format::MDRawBed, std::size_t > md_rawbeds_map;
        std::vector< ago::format::MDRawBed > md_rawbeds;
        ago::format::MDRawBed md_rawbed;

        // std::vector< AnnotationRawBed<> > annotation_rawbeds;
        // Sam2RawBed< std::vector< ago::format::MDSam<> >* > sam2bed;

        std::vector< std::pair<
            std::string,            // sam_lines
            ago::format::MDSam<>    // sams
        >> sam_pair_vec;

        for( size_t smp = 0; smp < fastq_paths.size(); ++smp )
        {
            load_counts = 0;

            sample_name = get_sample_name( fastq_paths[ smp ] );
            statistic_samples.emplace_back( sample_name, std::vector< double >
                    { 0.0, 0.0, 0.0, 0.0, 0.0 });
            //     RawReads FilteredReads Mappable%
            //          ReadsWithN  Mappable

            // print_mem_usage( "Sample " + sample_name + " Start" );

            monitor.set_monitor(( "	" + sample_name ).c_str(), 5 );
            monitor.log(( "	" + sample_name ).c_str(), "Aligning ..." );

            if( output_sam_ )
            {
                output_sam.open( db.output_dir().string() + sample_name + ".sam" );
                output_sam << sam_header_ss_.str();
            }

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

                    insert_statistic( rawread_count, fastq.getName() ); 

                    if( n_check( fastq ))
                    {
                        insert_statistic( fastq_n_count, fastq.getName() ); 
                        --job;
                        continue;
                    }

                    insert_statistic( fastq_count, fastq.getName() ); 

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

                parallel_pool.job_post([ fastq_str = fastq_ss.str(), &sam_pair_vec, &sam_mutex, load_counts, this ] () mutable
                {
                    aligning( fastq_str, sam_pair_vec, sam_mutex, load_counts );
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

            monitor.log(( "	" + sample_name ).c_str(), "Aligning ... Done" );
            monitor.log(( "	" + sample_name ).c_str(), "Converting ..." );

            // print_mem_usage( "Sample " + sample_name + " Aligning" );

            for( auto& sam_pair : sam_pair_vec )
            {
                if( output_sam_ )
                {
                    output_sam << sam_pair.first << "\n";
                }

                align_count.emplace( std::get< 0 >( sam_pair.second.data ));
                md_rawbed = ago::format::MDRawBed( sam_pair.second );
                md_rawbed.reducing_tail( max_tail_len );

                if( md_rawbeds_map.find( md_rawbed ) == md_rawbeds_map.end() )
                    md_rawbeds_map[ md_rawbed ] = 0;

                md_rawbeds_map[ md_rawbed ]++;

                // sams.emplace_back(   std::move( sam_pair.second ));
                // sam_pair.first = "";
            }

            sam_pair_vec.clear();
            output_sam.close();

            // auto raw_beds( sam2bed.Convert( &sams ));
            // sams.clear();

            // print_mem_usage( "Sample " + sample_name + " Sam Releasing" );

            // for( auto itr = raw_beds->begin(); itr != raw_beds->end(); ++itr )
            // {
            //     annotation_rawbeds.emplace_back( AnnotationRawBed<>( itr->first ));
            // }

            // print_mem_usage( "Sample " + sample_name + " Bed Converting" );

            // sam2bed.rawbed_map_->clear();
            // sam2bed.rawbed_map2_->clear();
            // raw_beds->clear();

            for( auto& mdbed : md_rawbeds_map )
            {
                (*(( uint32_t* )( mdbed.first.get_reads_count() ))) = mdbed.second;
                md_rawbeds.emplace_back( mdbed.first );
            }

            md_rawbeds_map.clear();
            db.bed_samples.emplace_back( sample_name, md_rawbeds );
            md_rawbeds.clear();

            // db.bed_samples.emplace_back( sample_name, annotation_rawbeds );
            // annotation_rawbeds.clear();

            // print_mem_usage( "Sample " + sample_name + " Bed Releasing" );
            
            for( auto& raw : rawread_count )
            {
                statistic_samples[ smp ].second[0] += raw.second;
            }
            
            for( auto& fqn : fastq_n_count )
            {
                statistic_samples[ smp ].second[1] += fqn.second;
            }

            for( auto& fqc : fastq_count )
            {
                statistic_samples[ smp ].second[2] += fqc.second;

                if( align_count.find( fqc.first ) != align_count.end() )
                {
                    statistic_samples[ smp ].second[3] += fqc.second;
                }
            }

            rawread_count.clear();
            fastq_n_count.clear();
            fastq_count.clear();
            align_count.clear();

            statistic_samples[ smp ].second[4] =
                statistic_samples[ smp ].second[3] * 100 / statistic_samples[ smp ].second[2];

            // print_mem_usage( "Sample " + sample_name + " Statistic Releasing" );

            monitor.log(( "	" + sample_name ).c_str(), "Converting ... Done" );
            monitor.log(( "	" + sample_name ).c_str(), "Done" );
            monitor.log( "Component GithubTailorFastqToBed", ( sample_name ).c_str() );

            // print_mem_usage( "Sample " + sample_name + " Complete" );
        }

        make_statistic( db.output_dir().string(), statistic_samples );
        monitor.log( "Component GithubTailorFastqToBed", "Complete" );

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

    void insert_statistic( std::map< std::string, size_t >& counts, const std::string& fq_name )
    {
        std::vector< std::string > split;
        boost::iter_split( split, fq_name, boost::algorithm::first_finder( " " ));

        if( counts.find( split[0] ) == counts.end() )
            counts[ split[0] ] = 0;
        
        counts[ split[0] ] += 1;
    }

    void make_statistic( auto& output_dir, auto& statistic_samples )
    {
        std::ofstream output( output_dir + "mappability.text" );

        output << "Sample";

        for( auto& sts : statistic_samples )
        {
            output << "\t" << sts.first;
        }

        output << "\n";

        for( size_t i = 0; i < 5; ++i )
        {
            switch( i )
            {
                case 0 : output << "RawReads:"; break;
                case 1 : output << "ReadsWithN:"; break;
                case 2 : output << "FilteredReads:"; break;
                case 3 : output << "Mappable:"; break;
                case 4 : output << "Mappable%:"; break;
                default: std::runtime_error( "out of row in mappability.text" );
            }

            for( auto& sts : statistic_samples )
            {
                output << std::fixed << std::setprecision( 0 ) << "\t" << sts.second[i];
            }
            
            output << "\n";
        }

        output.close();
    }

    void aligning( auto& fastq_str, auto& sam_pair_vec, auto& sam_mutex, auto& load_counts )
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
                abwtt_, &fastq_ss, &sam_ss, align_min_length_, align_min_multi_, align_allow_mismatch_, align_allow_t2c_, false 
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

        // SamDefaultTuple sam_tpl;
        // Sam<> sam;

        std::string samline;
        std::vector< std::string > split;

        std::vector< std::pair<
            std::string,            // sam_lines
            ago::format::MDSam<>    // sams
        >> sam_pairs;

        while( std::getline( sam_ss, samline ))
        {
            if( samline.at(0) == '@' ) continue;
            sam_pairs.emplace_back( std::make_pair( samline, ago::format::MDSam<>( samline, align_allow_t2c_ )));
        }

        // while( true )
        // {
        //     if( !std::getline( sam_ss, samline ))
        //     {
        //         break;
        //     }

        //     boost::iter_split( split, samline, boost::algorithm::first_finder( "\t" ));

        //     for( size_t x = 12; x < split.size(); ++x )
        //     {
        //         split[11] += ( '\t' + split[x] );
        //     }

        //     TupleUtility< SamDefaultTuple, std::tuple_size< SamDefaultTuple >::value -1 >::FillTuple( sam_tpl, split );
        //     std::get< 11 >( sam_tpl ) = typename std::tuple_element< 11, SamDefaultTuple >::type( split[ 11 ]); 

        //     sam.data = sam_tpl;
        //     sam_pairs.emplace_back( std::make_pair( samline, sam ));

        //     split.clear();
        // }

        // print_mem_usage( "Tailor Sam Converting-" + std::to_string( load_counts ));

        sam_ss.str("");
        sam_ss.clear();

        if( !sam_pairs.empty() )
        {
            std::lock_guard< std::mutex > sam_lock( sam_mutex );
            std::move( sam_pairs.begin(), sam_pairs.end(), std::back_inserter( sam_pair_vec ));
        }

        sam_pairs.clear();

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
