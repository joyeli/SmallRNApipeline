#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/analyzer_lendist.hpp>
#include <AGO/algorithm/analyzer_mirdist.hpp>

namespace ago {
namespace component {

void output_result_for_test( const std::string& name, const std::string& sample, const std::string& type,
        std::vector< std::map< std::string, std::map< std::string, double >>>& analyzer_result )
//         Each    LenMir_Map  Anno_Seed           All_Length      Length
{
    std::ofstream output( name + "_" + sample + "_" + type + ".test" );

    for( auto& map : analyzer_result )
    {
        for( auto& mir : map )
        {
            if( type != "" && mir.first != type && mir.first.at(0) != '.' )
                continue;

            for( auto& len : mir.second )
                output << sample << "\t" << mir.first << "\t" << len.first << "\t" << len.second << "\n";
        }
    }

    output.close();
}

class Analyzer : public engine::NamedComponent
{
    using Base = engine::NamedComponent;

  public:

    using Base::Base;

    virtual void initialize() override
    {
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        monitor.set_monitor( "Component Analyzer", db.bed_samples.size() +2);
        monitor.log( "Component Analyzer", "Start" );

        std::string output_path( db.output_dir().string() );

        ago::algorithm::AnalyzerLenDist len;
        ago::algorithm::AnalyzerMirDist mir;

        for( auto& sample : db.bed_samples )
        {
            len.analyzer.run( &sample.second, 0, true, output_path, sample.first, db.genome_table, db.analyzer_result, 0 );
            len.tailing_ratio( db.analyzer_result );


            mir.analyzer.run( &sample.second, 0, true, output_path, sample.first, db.genome_table, db.analyzer_result, 0 );
            mir.tailing_ratio( db.analyzer_result );
            mir.tailing_ratio_for_PM( db.analyzer_result, sample.second, db.genome_table );
            mir.tailing_ratio_for_miRNA( db.analyzer_result, sample.second, db.genome_table );
            mir.tailing_detail_for_miRNA( sample, db.genome_table, output_path );

            annotation_tailing_detail( sample, db.genome_table, output_path );
            un_annotation_tailing_detail( sample, db.genome_table, output_path );

            monitor.log( "Component Analyzer", ( sample.first ).c_str() );

            db.analyzer_result_samples.emplace_back( sample.first, db.analyzer_result );
            output_result_for_test( output_path + "AnaDone", sample.first, "miRNA", db.analyzer_result );
            db.analyzer_result.clear();
            sample.second.clear();
        }

        db.bed_samples.clear();

        monitor.log( "Component Analyzer", "Complete" );
    }

    void annotation_tailing_detail(
              std::pair< std::string, std::vector< AnnotationRawBed<> >>& sample_rawbed
            , std::map< std::string, std::string >& genome_table
            , std::string& output_path
    )
    {
        std::map< std::string, std::map< size_t, std::map< std::string, double >>> anno_tail_map;
        //        anno                    lens              tail_seq    count

        std::map< std::string, std::map< size_t, std::map< std::string, double >>>::iterator anno_tail_it;
        std::map< size_t, std::map< std::string, double >>::iterator len_tail_it;
        std::map< std::string, double >::iterator tail_it;

        for( auto& raw_bed : sample_rawbed.second )
        {
            for( auto& raw_bed_info : raw_bed.annotation_info_ )
            {
                for( int i = 0; i < raw_bed_info.size(); i+=2 )
                {
                    if( raw_bed_info[i] != "miRNA" )
                    {
                        std::string anno( raw_bed_info[i+1] + "|" + raw_bed_info[i] + "|" + std::to_string( raw_bed.is_filtered_ ));
                        std::string tail = raw_bed.getTail();
                        anno_cehck( anno_tail_map, raw_bed, anno, tail );
                    }
                }
            }
        }

        std::ofstream output( output_path + "/" + sample_rawbed.first + "_anno_tailing.tsv" );
        output << "Anno\tType\tFilter\tLen\tTailSeq\tCount\n";

        for( auto& anno_len : anno_tail_map )
        {
            for( auto& len_tail : anno_len.second )
            {
                for( auto& tail_count : len_tail.second )
                {
                    std::vector< std::string > anno_type;
                    boost::iter_split( anno_type, anno_len.first, boost::algorithm::first_finder( "|" ));

                    output
                        << anno_type[0] << "\t"
                        << anno_type[1] << "\t"
                        << anno_type[2] << "\t"
                        << len_tail.first << "\t"
                        << tail_count.first << "\t"
                        << tail_count.second << "\n"
                        ;
                }
            }
        }

        output.close();
    }

    void un_annotation_tailing_detail(
              std::pair< std::string, std::vector< AnnotationRawBed<> >>& sample_rawbed
            , std::map< std::string, std::string >& genome_table
            , std::string& output_path
    )
    {
        std::map< std::string, std::map< size_t, std::map< std::string, double >>> anno_tail_map;
        //        anno                    lens              tail_seq    count

        for( auto& raw_bed : sample_rawbed.second )
        {
            for( auto& raw_bed_info : raw_bed.annotation_info_ )
            {
                if( raw_bed_info.size() == 0 )
                {
                    std::string anno = raw_bed.chromosome_ + ":" + std::to_string( raw_bed.start_ + 1 ) + "-" + std::to_string( raw_bed.end_ );
                    std::string tail = raw_bed.getTail();
                    anno_cehck( anno_tail_map, raw_bed, anno, tail );
                }
            }
        }

        std::ofstream output( output_path + "/" + sample_rawbed.first + "_un_anno_tailing.tsv" );
        output << "Anno\tLen\tTailSeq\tCount\n";

        for( auto& anno_len : anno_tail_map )
        {
            for( auto& len_tail : anno_len.second )
            {
                for( auto& tail_count : len_tail.second )
                {
                    output
                        << anno_len.first << "\t"
                        << len_tail.first << "\t"
                        << tail_count.first << "\t"
                        << tail_count.second << "\n"
                        ;
                }
            }
        }

        output.close();
    }

    void anno_cehck( auto& anno_tail_map, const auto& raw_bed, const auto& anno, const auto& tail )
    {
        size_t length = raw_bed.length_ - raw_bed.tail_length_;
        double read_count = (double)(raw_bed.reads_count_) / (double)(raw_bed.multiple_alignment_site_count_);

        std::map< std::string, std::map< size_t, std::map< std::string, double >>>::iterator anno_tail_it = anno_tail_map.find( anno );

        if( anno_tail_it != anno_tail_map.end() )
        {
            std::map< size_t, std::map< std::string, double >>::iterator len_tail_it = anno_tail_it->second.find( length );

            if( len_tail_it != anno_tail_it->second.end() )
            {
                std::map< std::string, double >::iterator tail_it = len_tail_it->second.find( tail );

                if( tail_it != len_tail_it->second.end() )
                {
                    tail_it->second += read_count;
                }
                else
                {
                    len_tail_it->second.emplace( tail, read_count );
                }
            }
            else
            {
                std::map< std::string, double > tail_temp;
                tail_temp.emplace( tail, read_count );

                anno_tail_it->second.emplace( length, tail_temp );
            }
        }
        else
        {
            std::map< std::string, double > tail_temp;
            tail_temp.emplace( tail, read_count );

            std::map< size_t, std::map< std::string, double >> len_tail_temp;
            len_tail_temp.emplace( length, tail_temp );

            anno_tail_map.emplace( anno, len_tail_temp );
        }
    }
};

} // end of namespace component
} // end of namespace ago
