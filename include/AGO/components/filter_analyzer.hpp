#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/meta_analyzer.hpp>
#include <iomanip>

namespace ago {
namespace component {

class FilterAnalyzer : public engine::NamedComponent
{
    using DropTypeList = boost::mpl::vector<
        boost::mpl::vector< boost::mpl::string< 'rm', 'sk' >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>>;

    //                                   annotation           length, 0=sum  count
    using CountingTableType = std::map< std::string, std::map< std::size_t, double >>;
    using TailingTableType  = std::map< std::string, std::map< std::size_t, std::map< std::string, double >>>;
    //                                   annotation               length                  tail      count
    //                                       annotation_index             length_index
    using AnnoLengthIndexType = std::pair< std::set< std::string >, std::set< std::size_t >>;
    using BedSampleType = std::pair< std::string, std::vector< AnnotationRawBed<> >>;    
    using Filters = FilterWorker< AnnotationRawBed<>, DropTypeList >;
    using Base = engine::NamedComponent;

    std::vector< std::string > biotype_list;

    double sudo_count;
    bool is_filter;

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        sudo_count = p.get_optional< double >( "sudo_count" ).value_or( 0.000001 );
        is_filter  = p.get_optional< bool   >( "is_filter"  ).value_or( true     );

        if(  p.get_child_optional( "biotype_list" ))
        {
            for( auto& biotype : p.get_child( "biotype_list" ))
            {
                biotype_list.emplace_back( biotype.second.data() );
            }
        }
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

        monitor.set_monitor( "Component FilterAnalyzer", 4 + 11 * biotype_list.size() );
        monitor.log( "Component FilterAnalyzer", "Start" );
        auto bed_samples = db.bed_samples;

        monitor.log( "Component FilterAnalyzer", "Biotype Analysis ... " );
        output_biotype( db.output_dir().string(), db.genome_table, bed_samples );

        if( is_filter )
        {
            drop_filtering( bed_samples );
            output_biotype( db.output_dir().string(), db.genome_table, bed_samples, true );
        }

        AnnoLengthIndexType ano_len_idx;

        //                             gmpm          tailing_ratio
        std::vector< std::pair< CountingTableType, CountingTableType >> counting_tables;
        std::vector< std::pair< CountingTableType, CountingTableType >> ppm_counting_tables;

        for( std::size_t i = 0; i < biotype_list.size(); ++i )
        {
            auto& biotype = biotype_list[i];

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Get Index of Annotation and Length ... " );
            ano_len_idx = get_ano_len_idx( db.genome_table, bed_samples, biotype );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Get Counting Tables ... " );
            counting_tables = get_counting_tables( db.genome_table, bed_samples, ano_len_idx, biotype );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] PPM Normalization ... " );
            ppm_counting_tables = ppm_counting_tables_converter( counting_tables );

            boost::filesystem::create_directory( boost::filesystem::path( db.output_dir().string() + "/" + biotype ));

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Counting Tables ... " );
            output_biotype_detail( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, counting_tables, "count" );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing PPM Tables ... " );
            output_biotype_detail( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables, "ppm" );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Quantile Normalization ... " );
            quantile_normalize( ano_len_idx, ppm_counting_tables );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Quantiled PPM ... " );
            output_biotype_detail( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables, "quantile" );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Loading Difference ... " );
            output_loading_difference( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Length Difference ... " );
            output_length_difference( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Arms Difference ... " );
            if( biotype == "miRNA" )
                output_arms_difference( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables );

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Annotation Tailing ... " );
            output_annotated_tailing( db.output_dir().string() + "/" + biotype, db.genome_table, bed_samples, biotype );
        }

        output_biotype( db.output_dir().string(), db.genome_table, bed_samples, true, true );

        monitor.log( "Component FilterAnalyzer", "Outputing Non Annotation Tailing ... " );
        output_non_annotated_tailing( db.output_dir().string(), bed_samples );

        monitor.log( "Component FilterAnalyzer", "Complete" );
    }

    void drop_filtering( std::vector< BedSampleType >& bed_samples )
    {
		Filters run_filter;
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &run_filter, &bed_samples ] () mutable
            {
                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    anno_rawbed = run_filter.Filter( anno_rawbed );
                }
            });
        }
        smp_parallel_pool.flush_pool();
    }

    AnnoLengthIndexType get_ano_len_idx(
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            const std::string& biotype = "" // default annotation type
            )
    {
        std::vector< std::set< std::string >> smp_ano_idx( bed_samples.size(), std::set< std::string >() );
        std::vector< std::set< std::size_t >> smp_len_idx( bed_samples.size(), std::set< std::size_t >() );

        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &genome_table, &biotype, &bed_samples, &smp_ano_idx, &smp_len_idx, this ] () mutable
            {
                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    if( raw_bed.is_filtered_ == 1 )
                        continue;

                    smp_len_idx[ smp ].emplace( raw_bed.length_ - raw_bed.tail_length_ );

                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                        {
                            if( biotype != "" && raw_bed_info[i] != biotype ) continue;
                            smp_ano_idx[ smp ].emplace( raw_bed_info[ i + ( biotype == "" ? 0 : 1 )]
                                + ( biotype == "" ? "" : ( "_" + raw_bed.getReadSeq( genome_table ).substr( 1, 7 )) ));
                        }
                    }
                }
            });
        }

        smp_parallel_pool.flush_pool();

        std::set< std::string > ano_idx;
        std::set< std::size_t > len_idx;
        std::set< std::size_t > len_temp;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            ano_idx.insert( smp_ano_idx[ smp ].begin(), smp_ano_idx[ smp ].end() );
            len_idx.insert( smp_len_idx[ smp ].begin(), smp_len_idx[ smp ].end() );
        }

        std::size_t count = 0;
        std::size_t last_len;

        for( auto& len : len_idx )
        {
            if( count == 0 ) last_len = len -1;
            if( len - last_len > 1 )
                for( std::size_t i = last_len + 1; i < len; ++i )
                    len_temp.emplace( i ); 
            count++;
        }

        len_idx.insert( len_temp.begin(), len_temp.end() );

        return std::make_pair( ano_idx, len_idx );
    }

    void output_biotype(
            const std::string& output_path, 
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            const bool& is_dropped = false,
            const bool& is_quantile = false
    )
    {
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        std::vector< CountingTableType > anno_table( bed_samples.size(), CountingTableType() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &bed_samples, &anno_table ] () mutable
            {
                //          annotation  count
                std::map< std::string, double > anno_check;
                double anno_check_sum;
                std::size_t read_len;

                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    if( raw_bed.is_filtered_ == 1 )
                        continue;

                    anno_check.clear();
                    anno_check_sum = 0.0;
                    read_len = raw_bed.length_ - raw_bed.tail_length_;

                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                        {
                            anno_check_sum++;
                            anno_check[ raw_bed_info[i] ]
                                = anno_check.find( raw_bed_info[i] ) != anno_check.end()
                                ? anno_check[ raw_bed_info[i] ] + 1.0 : 1.0;
                        }
                    }

                    for( auto& anno : anno_check )
                    {
                        anno.second = anno.second / anno_check_sum;

                        if( anno_table[ smp ].find( anno.first ) != anno_table[ smp ].end() )
                        {
                            anno_table[ smp ][ anno.first ][ read_len ]
                                = anno_table[ smp ][ anno.first ].find( read_len ) != anno_table[ smp ][ anno.first ].end()
                                ? anno_table[ smp ][ anno.first ][ read_len ] + anno.second
                                : anno.second
                                ;
                        }
                        else
                        {
                            anno_table[ smp ][ anno.first ] = std::map< std::size_t, double >();
                            anno_table[ smp ][ anno.first ][ read_len ] = anno.second;
                        }
                    }
                }
            });
        }

        smp_parallel_pool.flush_pool();
        AnnoLengthIndexType ano_len_idx = get_ano_len_idx( genome_table, bed_samples );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &output_path, &is_dropped, &is_quantile, sample_name = bed_samples[ smp ].first, &ano_len_idx, &anno_table ] () mutable
            {
                std::ofstream output( output_path + "/" + sample_name + "_biotype"
                    + ( is_dropped  ? "_filtered"  : "" )
                    + ( is_quantile ? "_quantiled" : "" )
                    + ".tsv" );

                output << "Annotation\tSum";

                std::string out_temp = "";
                double sum = 0.0;

                for( auto& len : ano_len_idx.second )
                {
                    output << "\t" << len;
                }

                for( auto& anno : ano_len_idx.first )
                {
                    output << "\n" << anno;

                    out_temp = "";
                    sum = 0.0;

                    if( anno_table[ smp ].find( anno ) != anno_table[ smp ].end() )
                    {
                        for( auto& len : ano_len_idx.second )
                        {
                            if( anno_table[ smp ][ anno ].find( len ) != anno_table[ smp ][ anno ].end() )
                            {
                                out_temp += "\t" + std::to_string( anno_table[ smp ][ anno ][ len ] );
                                sum += anno_table[ smp ][ anno ][ len ];
                            }
                            else
                            {
                                out_temp += "\t0";
                            }
                        }
                    }
                    else
                    {
                        for( auto& len : ano_len_idx.second )
                            out_temp += "\t0";
                    }

                    output << "\t" << sum << out_temp;
                }

                output << "\n";
                output.close();
            });
        }

        smp_parallel_pool.flush_pool();
    }

    std::vector< std::pair< CountingTableType, CountingTableType >> get_counting_tables(
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            const std::string& biotype
            )
    {
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        std::vector< std::pair< CountingTableType, CountingTableType >> counting_tables(
                bed_samples.size(), std::pair< CountingTableType, CountingTableType >()
                );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &genome_table, &bed_samples, &biotype, &ano_len_idx, &counting_tables, this ] () mutable
            {
                CountingTableType total_table;
                CountingTableType tail_table;

                std::size_t tail_len = 0;
                std::size_t length = 0;

                double count = 0.0;
                std::set< std::string > annos;

                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    if( raw_bed.is_filtered_ == 1 )
                        continue;

                    annos.clear();

                    tail_len = raw_bed.tail_length_; 
                    length   = raw_bed.length_ - raw_bed.tail_length_;
                    count    = (double)(raw_bed.reads_count_) / (double)(raw_bed.multiple_alignment_site_count_);

                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                        {
                            if( raw_bed_info[i] != biotype ) continue;
                            annos.emplace( raw_bed_info[ i+1 ] + "_" + raw_bed.getReadSeq( genome_table ).substr( 1, 7 ));
                        }
                    }

                    for( auto& annotation : annos )
                    {
                        counting_table_inserter( total_table, annotation, length, count + sudo_count );

                        if( tail_len != 0 )
                            counting_table_inserter( tail_table, annotation, length, count + sudo_count );
                    }
                }
                
                counting_table_refinder( total_table, ano_len_idx, sudo_count );
                counting_table_refinder( tail_table, ano_len_idx, sudo_count );
                tailing_ratio_transfer( total_table, tail_table, sudo_count );

                counting_tables[ smp ] = std::make_pair( total_table, tail_table );
            });
        }
        smp_parallel_pool.flush_pool();
        return counting_tables;
    }

    std::vector< std::pair< CountingTableType, CountingTableType >> ppm_counting_tables_converter(
            const std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        ParaThreadPool smp_parallel_pool( counting_tables.size() );
        std::vector< std::pair< CountingTableType, CountingTableType >> ppm_counting_tables(
                counting_tables.size(), std::pair< CountingTableType, CountingTableType >()
                );

        for( std::size_t smp = 0; smp < counting_tables.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &ppm_counting_tables, &counting_tables ] () mutable
            {
                CountingTableType ppm_table = counting_tables[ smp ].first;

                double seq_sum = 0.0;
                double ppm_ratio = 1000000;

                for( auto& anno : ppm_table )
                {
                    seq_sum += anno.second[0];
                }

                ppm_ratio = ppm_ratio / seq_sum;

                for( auto& anno : ppm_table )
                {
                    for( auto& len : anno.second )
                    {
                        len.second = len.second * ppm_ratio;
                    }
                }

                ppm_counting_tables[ smp ] = std::make_pair( ppm_table, std::move( counting_tables[ smp ].second ));
            });
        }
        smp_parallel_pool.flush_pool();
        return ppm_counting_tables;
    }

    void counting_table_inserter(
            CountingTableType& counting_table,
            const std::string& anno,
            const std::size_t& length,
            const double& count
            )
    {
        if( counting_table.find( anno ) != counting_table.end() )
        {
            if( counting_table[ anno ].find( length ) != counting_table[ anno ].end() )
            {
                counting_table[ anno ][ length ] += count;
            }
            else
            {
                counting_table[ anno ][ length ] = count;
            }
        }
        else
        {
            std::map< std::size_t, double > temp;
            temp[ length ] = count;
            counting_table[ anno ] = temp;
        }
    }

    void counting_table_refinder( CountingTableType& counting_table, const AnnoLengthIndexType& ano_len_idx, const double& sudo_count = 0.0 )
    {
        for( auto& anno : ano_len_idx.first )
        {
            if( counting_table.find( anno ) != counting_table.end() )
            {
                for( auto& len : ano_len_idx.second )
                {
                    if( counting_table[ anno ].find( len ) == counting_table[ anno ].end() )
                    {
                        counting_table[ anno ][ len ] = sudo_count;
                    }
                }
            }
            else
            {
                std::map< std::size_t, double > temp;
                for( auto& len : ano_len_idx.second )
                {
                    temp[ len ] = sudo_count;
                }
                counting_table[ anno ] = temp;
            }
        }

        double sum = 0.0;

        for( auto& anno : counting_table )
        {
            sum = 0.0;

            for( auto& len : anno.second )
            {
                sum += len.second;
            }

            anno.second[0] = sum;
        }
    }

    void tailing_ratio_transfer( CountingTableType& total_table, CountingTableType& tail_table, const double& sudo_count = 0.0 )
    {
        for( auto& anno : tail_table )
        {
            for( auto& len : anno.second )
            {
                if( len.second != 0 && len.second != sudo_count )
                {
                    len.second = len.second / total_table[ anno.first ][ len.first ];
                }
            }
        }
    }

    void quantile_normalize(
            AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::vector< QuantileDataType<> > ppm_qvec;
        std::vector< double > ppm_vec;

        for( std::size_t smp = 0; smp < counting_tables.size(); ++smp )
        {
            auto& counting_table = counting_tables[ smp ].first;
            ppm_vec.clear();

            for( auto& anno : ano_len_idx.first )
            {
                ppm_vec.emplace_back( counting_table[ anno ][0] );

                for( auto& len : ano_len_idx.second )
                {
                    counting_table[ anno ][ len ]
                        = counting_table[ anno ][ len ]
                        / counting_table[ anno ][0]
                        ;
                }
            }

            ppm_qvec.emplace_back( QuantileDataType<>( ppm_vec ));
        }

        QuantileNor quntile( ppm_qvec );
        std::size_t qidx = 0;

        for( std::size_t smp = 0; smp < counting_tables.size(); ++smp )
        {
            auto& counting_table = counting_tables[ smp ].first;
            qidx = 0;

            for( auto& anno : ano_len_idx.first )
            {
                counting_table[ anno ][0] = ppm_qvec[ smp ].value_[ qidx ].first;

                for( auto& len : ano_len_idx.second )
                {
                    counting_table[ anno ][ len ]
                        = counting_table[ anno ][ len ]
                        * counting_table[ anno ][0]
                        ;
                }
                qidx++;
            }
        }
    }

    void output_biotype_detail(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables,
            const std::string& tag
            )
    {
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &output_path, &bed_samples, &ano_len_idx, &counting_tables, &tag ] () mutable
            {
                double pm = 0.0;

                std::ofstream out_gmpm( output_path + "/" + bed_samples[ smp ].first + "_" + tag + "_GMPM.tsv" );
                std::ofstream out_gm  ( output_path + "/" + bed_samples[ smp ].first + "_" + tag + "_GM.tsv" );
                std::ofstream out_pm  ( output_path + "/" + bed_samples[ smp ].first + "_" + tag + "_PM.tsv" );
                std::ofstream out_tail( output_path + "/" + bed_samples[ smp ].first + "_" + tag + "_TailRatio.tsv" );

                out_gmpm << "Annotation" ;// \tSum";
                out_gm   << "Annotation" ;// \tSum";
                out_pm   << "Annotation" ;// \tSum";
                out_tail << "Annotation" ;// \tSum";

                for( auto& len : ano_len_idx.second )
                {
                    out_gmpm << "\t" << len;
                    out_gm   << "\t" << len;
                    out_pm   << "\t" << len;
                    out_tail << "\t" << len;
                }

                for( auto& anno : ano_len_idx.first )
                {
                    pm = counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];

                    out_gmpm << "\n" << anno ;// << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].first[ anno ][ 0 ];
                    out_gm   << "\n" << anno ;// << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].first[ anno ][ 0 ] - pm;
                    out_pm   << "\n" << anno ;// << "\t" << std::fixed << std::setprecision( 2 ) << pm;
                    out_tail << "\n" << anno ;// << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].second[ anno ][ 0 ];

                    for( auto& len : ano_len_idx.second )
                    {
                        pm = counting_tables[ smp ].first[ anno ][ len ] * counting_tables[ smp ].second[ anno ][ len ];

                        out_gmpm << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].first[ anno ][ len ];
                        out_gm   << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].first[ anno ][ len ] - pm;
                        out_pm   << "\t" << std::fixed << std::setprecision( 2 ) << pm;
                        out_tail << "\t" << std::fixed << std::setprecision( 2 ) << counting_tables[ smp ].second[ anno ][ len ];
                    }
                }

                out_gmpm << "\n";
                out_gm   << "\n";
                out_pm   << "\n";
                out_tail << "\n";

                out_gmpm.close();
                out_gm.close();
                out_pm.close();
                out_tail.close();
            });
        }
        smp_parallel_pool.flush_pool();
    }

    double get_sum( const std::vector< double >& vec )
    {
        double res = 0;
        for( auto& val : vec )
            res += val;
        return res;
    }

    std::tuple< double, std::size_t, std::size_t > get_difference( const std::vector< double >& vec )
    {
        double sum = get_sum( vec );
        std::pair< double, std::size_t > max;
        std::pair< double, std::size_t > min;


        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == 0 )
            {
                max = std::make_pair( vec[ smp ], smp );
                min = std::make_pair( vec[ smp ], smp );
                continue;
            }

            if( max.first < vec[ smp ] )
            {
                max = std::make_pair( vec[ smp ], smp );
            }

            if( min.first > vec[ smp ] )
            {
                min = std::make_pair( vec[ smp ], smp );
            }
        }

        std::vector< double > folds;
        double fold = (( max.first/sum )-( min.first/sum ))/( min.first/sum );

        for( std::size_t smp = 0; smp < vec.size(); ++smp )
        {
            if( smp == max.second || smp == min.second ) continue;
            folds.emplace_back( (( vec[ smp ]/sum )-( min.first/sum ))/( min.first/sum ));
        }

        bool fold_tag = false;

        for( auto& fd : folds )
            if( fd >= 1 ) fold_tag = true;

        return{ fold_tag ? ( -1 * fold ) : fold, max.second, min.second };
    }

    std::vector< std::size_t > get_lengths( const std::vector< std::map< std::size_t, double >>& len_vec )
    {
        double max;
        std::size_t len;
        std::vector< std::size_t > res( len_vec.size(), std::size_t() );

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            len = 0;
            max = 0.0;

            for( auto& lens : len_vec[ smp ] )
            {
                if( max < lens.second )
                {
                    len = lens.first;
                    max = lens.second;
                }
            }

            res[ smp ] = len;
        }

        return res;
    }

    std::tuple< std::size_t, std::size_t, std::size_t > get_length_difference( const std::vector< std::size_t >& len_vec )
    {
        std::pair< std::size_t, std::size_t > max;
        std::pair< std::size_t, std::size_t > min;

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            if( smp == 0 )
            {
                max = std::make_pair( len_vec[ smp ], smp );
                min = std::make_pair( len_vec[ smp ], smp );
                continue;
            }

            if( max.first < len_vec[ smp ] )
            {
                max = std::make_pair( len_vec[ smp ], smp );
            }

            if( min.first > len_vec[ smp ] )
            {
                min = std::make_pair( len_vec[ smp ], smp );
            }
        }

        return { max.first - min.first, max.second, min.second };
    }

    std::vector< std::pair< std::string, double >> get_arms( const std::vector< std::map< std::string, double >>& arm_vec )
    {
        std::vector< std::pair< std::string, double >> res;
        std::pair< std::string, double > tmp = std::make_pair( "", 0.0 );

        for( std::size_t smp = 0; smp < arm_vec.size(); ++smp )
        {
            for( auto& arm : arm_vec[ smp ] )
            {
                if( tmp.second == 0.0 ) tmp = arm;
                if( tmp.second < arm.second ) tmp = arm;
            }

            res.emplace_back( tmp );
            tmp = std::make_pair( "", 0.0 );
        }

        return res;
    }

    std::string get_arm_difference( const std::vector< std::pair< std::string, double >>& arm_vec )
    {
        std::string res = "Y";
        std::set< std::string > samp;
        for( auto& arm : arm_vec ) samp.emplace( arm.first );
        if( samp.size() == 1 ) res = "N";
        return res;
    }

    void sort_difference( std::vector< std::pair< double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
            { return a.first > b.first; });
    }

    void output_loading_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_loading( output_path + "/LoadingDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/LoadingDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";
        out_tailing << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2";

        std::vector< double > vec_loading;
        std::vector< double > vec_tailing;

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::tuple< double, std::size_t, std::size_t > df_tuple_loading;
        std::tuple< double, std::size_t, std::size_t > df_tuple_tailing;

        double pm = 0.0;
        double gmpm = 0.0;
        double sum_loading = 0.0;

        std::string res_loading;
        std::string res_tailing;

        for( auto& anno : ano_len_idx.first )
        {
            vec_loading.clear();
            vec_tailing.clear();

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                pm   = counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];
                gmpm = counting_tables[ smp ].first[ anno ][ 0 ];

                vec_loading.emplace_back( gmpm );
                vec_tailing.emplace_back(( pm < 1 ? 0 : pm ) / gmpm * 100 );
            }

            sum_loading = get_sum( vec_loading );

            df_tuple_loading = get_difference( vec_loading );
            df_tuple_tailing = get_difference( vec_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_loading += std::to_string( std::get<0>( df_tuple_loading )) + "\t" + bed_samples[ std::get<1>( df_tuple_loading )].first + ":" + bed_samples[ std::get<2>( df_tuple_loading )].first;
            res_tailing += std::to_string( std::get<0>( df_tuple_tailing )) + "\t" + bed_samples[ std::get<1>( df_tuple_tailing )].first + ":" + bed_samples[ std::get<2>( df_tuple_tailing )].first;

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }

    void output_length_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_loading( output_path + "/LengthDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/LengthDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";
        out_tailing << "Annotation\tTotal\tLengthDifference\tSample1:Sample2\tLength1:Length2";

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::vector< std::size_t > len_loading;
        std::vector< std::size_t > len_tailing;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::size_t, std::tuple< double, double >>>> smp_anno_len_gmpms;

        std::vector< std::map< std::size_t, double >> vec_len_loading( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_tailing( bed_samples.size(), std::map< std::size_t, double >() );

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_loading;
        std::tuple< std::size_t, std::size_t, std::size_t > df_tuple_tailing;

        std::string res_loading;
        std::string res_tailing;

        std::vector< std::string > anno_splits; 
        std::string anno_tag;

        for( auto& anno : ano_len_idx.first )
        {
            anno_splits.clear();

            boost::iter_split( anno_splits, anno, boost::algorithm::first_finder( "_" ));
            anno_tag = anno_splits[0];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "_" + anno_splits[i] );
            anno_index.emplace( anno_tag );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( smp_anno_len_gmpms[ smp ][ anno_tag ].find( 0 ) == smp_anno_len_gmpms[ smp ][ anno_tag ].end() )
                    smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] = { 0, 0 };

                /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ 0 ] ) += counting_tables[ smp ].first[ anno ][ 0 ];

                for( auto& len : ano_len_idx.second )
                {
                    if( smp_anno_len_gmpms[ smp ][ anno_tag ].find( len ) == smp_anno_len_gmpms[ smp ][ anno_tag ].end() )
                        smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] = { 0, 0 };

                    /*GMPM*/ std::get<0>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) +=   counting_tables[ smp ].first[ anno ][ len ];
                    /* PM */ std::get<1>( smp_anno_len_gmpms[ smp ][ anno_tag ][ len ] ) += ( counting_tables[ smp ].first[ anno ][ len ] * counting_tables[ smp ].second[ anno ][ len ] );
                }
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_len_loading = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );
            vec_len_tailing = std::vector< std::map< std::size_t, double >>( bed_samples.size(), std::map< std::size_t, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& len : smp_anno_len_gmpms[ smp ][ anno ] )
                {
                    if( len.first == 0 )
                    {
                        vec_loading.emplace_back( std::get<0>( len.second ) );
                        continue;
                    }

                    vec_len_loading[ smp ][ len.first ] = std::get<0>( len.second );
                    vec_len_tailing[ smp ][ len.first ] = std::get<1>( len.second );
                }
            }

            sum_loading = get_sum( vec_loading );

            len_loading = get_lengths( vec_len_loading );
            len_tailing = get_lengths( vec_len_tailing );

            df_tuple_loading = get_length_difference( len_loading );
            df_tuple_tailing = get_length_difference( len_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t";

            res_loading += std::to_string( std::get<0>( df_tuple_loading )) + "\t";
            res_tailing += std::to_string( std::get<0>( df_tuple_tailing )) + "\t";

            res_loading += bed_samples[ std::get<1>( df_tuple_loading )].first + ":" + bed_samples[ std::get<2>( df_tuple_loading )].first + "\t";
            res_tailing += bed_samples[ std::get<1>( df_tuple_tailing )].first + ":" + bed_samples[ std::get<2>( df_tuple_tailing )].first + "\t";

            res_loading += std::to_string( len_loading[ std::get<1>( df_tuple_loading )]) + ":" + std::to_string( len_loading[ std::get<2>( df_tuple_loading )]);
            res_tailing += std::to_string( len_tailing[ std::get<1>( df_tuple_tailing )]) + ":" + std::to_string( len_tailing[ std::get<2>( df_tuple_tailing )]);

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }

    void output_arms_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_loading( output_path + "/ArmDifference_Loading.tsv" );
        std::ofstream out_tailing( output_path + "/ArmDifference_Tailing.tsv"   );

        out_loading << "Annotation\tTotal\tisLoadingDifference";
        out_tailing << "Annotation\tTotal\tisLoadingDifference";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_loading << "\tSample1";
            else out_loading << ":Sample" << smp + 1 ;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            if( smp == 0 ) out_tailing << "\tArm1";
            else out_tailing << ":Arm" << smp + 1 ;

        double sum_loading = 0.0;
        std::vector< double > vec_loading;

        std::set< std::string > anno_index;
        std::map< std::size_t, std::map< std::string, std::map< std::string, std::tuple< double, double >>>> smp_anno_arm_gmpms;

        std::vector< std::map< std::string, double >> vec_arm_loading( bed_samples.size(), std::map< std::string, double >() );
        std::vector< std::map< std::string, double >> vec_arm_tailing( bed_samples.size(), std::map< std::string, double >() );

        std::vector< std::pair< std::string, double >> arm_loading;
        std::vector< std::pair< std::string, double >> arm_tailing;

        std::vector< std::pair< double, std::string >> temp_loading;
        std::vector< std::pair< double, std::string >> temp_tailing;

        std::string is_diff_loading;
        std::string is_diff_tailing;

        std::string res_loading;
        std::string res_tailing;

        std::vector< std::string > anno_splits; 
        std::string anno_tag;
        std::string anno_arm;

        for( auto& anno : ano_len_idx.first )
        {
            anno_splits.clear();

            boost::iter_split( anno_splits, anno, boost::algorithm::first_finder( "_" ));
            anno_tag = anno_splits[0];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "_" + anno_splits[i] );
            anno_splits.clear();

            boost::iter_split( anno_splits, anno_tag, boost::algorithm::first_finder( "-" ));

            anno_tag = anno_splits[0];
            anno_arm = anno_splits[ anno_splits.size() -1 ];

            for( std::size_t i = 1; i < anno_splits.size() -1; ++i ) anno_tag += ( "-" + anno_splits[i] );
            anno_index.emplace( anno_tag );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                if( smp_anno_arm_gmpms[ smp ][ anno_tag ].find( anno_arm ) == smp_anno_arm_gmpms[ smp ][ anno_tag ].end() )
                    smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] = { 0, 0 };

                /*GMPM*/ std::get<0>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += counting_tables[ smp ].first[ anno ][ 0 ];
                /* PM */ std::get<1>( smp_anno_arm_gmpms[ smp ][ anno_tag ][ anno_arm ] ) += counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];
            }
        }

        for( auto& anno : anno_index )
        {
            vec_loading.clear();
            vec_arm_loading = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );
            vec_arm_tailing = std::vector< std::map< std::string, double >>( bed_samples.size(), std::map< std::string, double >() );

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                for( auto& arm : smp_anno_arm_gmpms[ smp ][ anno ] )
                {
                    vec_loading.emplace_back( std::get<0>( arm.second ) );
                    vec_arm_loading[ smp ][ arm.first ] =   std::get<0>( arm.second );
                    vec_arm_tailing[ smp ][ arm.first ] = ( std::get<1>( arm.second ) < 1 ? 0 : std::get<1>( arm.second )) / std::get<0>( arm.second ) * 100;
                }
            }

            sum_loading = get_sum( vec_loading );

            arm_loading = get_arms( vec_arm_loading );
            arm_tailing = get_arms( vec_arm_tailing );

            is_diff_loading = get_arm_difference( arm_loading );
            is_diff_tailing = get_arm_difference( arm_tailing );

            res_loading = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_loading;
            res_tailing = "\n" + anno + "\t" + std::to_string( sum_loading ) + "\t" + is_diff_tailing;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_loading += "\t" + bed_samples[ smp ].first;
                else res_loading += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_loading += "\t" + arm_loading[ smp ].first;
                else res_loading += ":" + arm_loading[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tailing += "\t" + bed_samples[ smp ].first;
                else res_tailing += ":" + bed_samples[ smp ].first;

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                if( smp == 0 ) res_tailing += "\t" + arm_tailing[ smp ].first;
                else res_tailing += ":" + arm_tailing[ smp ].first;

            temp_loading.emplace_back( std::make_pair( sum_loading, res_loading ));
            temp_tailing.emplace_back( std::make_pair( sum_loading, res_tailing ));
        }

        sort_difference( temp_loading );
        sort_difference( temp_tailing );

        for( auto& output : temp_loading ) out_loading << output.second;
        for( auto& output : temp_tailing ) out_tailing << output.second;

        out_loading << "\n";
        out_tailing << "\n";

        out_loading.close();
        out_tailing.close();
    }

    void output_annotated_tailing(
            const std::string& output_path, 
            std::map< std::string, std::string >& genome_table,
            std::vector< BedSampleType >& bed_samples,
            const std::string& biotype
            )
    {
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &genome_table, &biotype, &output_path, &bed_samples, this ] () mutable
            {
                TailingTableType anno_tail_map;

                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    if( raw_bed.is_filtered_ == 1 )
                        continue;

                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                        {
                            if( raw_bed_info[i] != biotype )
                                continue;

                            anno_tail_inserter(
                                anno_tail_map,
                                raw_bed_info[ i+1 ] + "_" + raw_bed.getReadSeq( genome_table ).substr( 1, 7 ),
                                raw_bed.getTail(),
                                raw_bed.length_ - raw_bed.tail_length_,
                                (double)(raw_bed.reads_count_) / (double)(raw_bed.multiple_alignment_site_count_)
                                );
                        }
                    }
                }
                output_tailing_table( output_path + "/" + bed_samples[ smp ].first, anno_tail_map );
            });
        }
        smp_parallel_pool.flush_pool();
    }

    void output_non_annotated_tailing(
            const std::string& output_path, 
            std::vector< BedSampleType >& bed_samples
            )
    {
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &output_path, &bed_samples, this ] () mutable
            {
                TailingTableType anno_tail_map;

                for( auto& raw_bed : bed_samples[ smp ].second )
                {
                    for( auto& raw_bed_info : raw_bed.annotation_info_ )
                    {
                        if( raw_bed_info.size() != 0 )
                            continue;

                        anno_tail_inserter(
                            anno_tail_map,
                            raw_bed.chromosome_ + ":" + std::to_string( raw_bed.start_ + 1 ) + "-" + std::to_string( raw_bed.end_ ),
                            raw_bed.getTail(),
                            raw_bed.length_ - raw_bed.tail_length_,
                            (double)(raw_bed.reads_count_) / (double)(raw_bed.multiple_alignment_site_count_)
                            );
                    }
                }
                output_tailing_table( output_path + "/" + bed_samples[ smp ].first + "_non-annotated", anno_tail_map );
            });
        }
        smp_parallel_pool.flush_pool();
    }

    void anno_tail_inserter(
            TailingTableType& anno_tail_map,
            const std::string& anno,
            const std::string& tail,
            const std::size_t& length,
            const double& read_count
            )
    {
        if( anno_tail_map.find( anno ) != anno_tail_map.end() )
        {
            if( anno_tail_map[ anno ].find( length ) != anno_tail_map[ anno ].end() )
            {
                if( anno_tail_map[ anno ][ length ].find( tail ) != anno_tail_map[ anno ][ length ].end() )
                {
                    anno_tail_map[ anno ][ length ][ tail ] += read_count;
                }
                else
                {
                    anno_tail_map[ anno ][ length ][ tail ] = read_count;
                }
            }
            else
            {
                std::map< std::string, double > tail_temp;

                tail_temp[ tail ] = read_count;
                anno_tail_map[ anno ][ length ] = tail_temp;
            }
        }
        else
        {
            std::map< std::string, double > tail_temp;
            std::map< size_t, std::map< std::string, double >> len_tail_temp;

            tail_temp[ tail ] = read_count;
            len_tail_temp[ length ] = tail_temp;
            anno_tail_map[ anno ] = len_tail_temp;
        }
    }

    void output_tailing_table(
            const std::string& output_path,
            const TailingTableType& anno_tail_map
            )
    {
        std::ofstream output( output_path + "_tailing.tsv" );
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
};

} // end of namespace component
} // end of namespace ago
