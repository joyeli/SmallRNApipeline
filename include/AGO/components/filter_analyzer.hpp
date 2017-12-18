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

        monitor.set_monitor( "Component FilterAnalyzer", 4 + 9 * biotype_list.size() );
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

            monitor.log( "Component FilterAnalyzer", "[ " + std::to_string( i ) + " / " + std::to_string( biotype_list.size() ) + " ][ " + biotype + " ] Outputing Sample Difference ... " );
            output_sample_difference( db.output_dir().string() + "/" + biotype, bed_samples, ano_len_idx, ppm_counting_tables );

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
                            counting_table_inserter( tail_table, annotation, length, count );
                    }
                }
                
                counting_table_refinder( total_table, ano_len_idx, sudo_count );
                counting_table_refinder( tail_table, ano_len_idx, 0.0 );
                tailing_ratio_transfer( total_table, tail_table );

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

    void tailing_ratio_transfer( CountingTableType& total_table, CountingTableType& tail_table )
    {
        for( auto& anno : tail_table )
        {
            for( auto& len : anno.second )
            {
                if( len.second != 0 )
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

        return{(( max.first/sum )-( min.first/sum ))/( min.first/sum ), max.second, min.second };
    }

    std::vector< std::size_t > get_length_difference( const std::vector< std::map< std::size_t, double >>& len_vec )
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

    std::size_t get_length_distance( const std::vector< std::size_t >& len_vec )
    {
        std::size_t max;
        std::size_t min;

        for( std::size_t smp = 0; smp < len_vec.size(); ++smp )
        {
            if( smp == 0 )
            {
                max = len_vec[ smp ];
                min = len_vec[ smp ];
                continue;
            }

            if( max < len_vec[ smp ] )
            {
                max = len_vec[ smp ];
            }

            if( min > len_vec[ smp ] )
            {
                min = len_vec[ smp ];
            }
        }

        return max - min;
    }

    void sort_difference( std::vector< std::pair< double, std::string >>& vec )
    {
        std::sort( vec.begin(), vec.end(),
            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
            { return a.first > b.first; });
    }

    void output_sample_difference(
            const std::string& output_path,
            const std::vector< BedSampleType >& bed_samples,
            const AnnoLengthIndexType& ano_len_idx,
            std::vector< std::pair< CountingTableType, CountingTableType >>& counting_tables
            )
    {
        std::ofstream out_gmpm( output_path + "/SampleDifference_GMPM.tsv" );
        std::ofstream out_gm  ( output_path + "/SampleDifference_GM.tsv"   );
        std::ofstream out_pm  ( output_path + "/SampleDifference_PM.tsv"   );

        out_gmpm << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2\tLengthDifference";
        out_gm   << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2\tLengthDifference";
        out_pm   << "Annotation\tTotal\tExpressionDifference\tSample1:Sample2\tLengthDifference";

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            out_gmpm << "\t" << bed_samples[ smp ].first << "Length";
            out_gm   << "\t" << bed_samples[ smp ].first << "Length";
            out_pm   << "\t" << bed_samples[ smp ].first << "Length";
        }

        std::vector< double > vec_gmpm;
        std::vector< double > vec_gm;
        std::vector< double > vec_pm;

        std::vector< std::size_t > len_gmpm;
        std::vector< std::size_t > len_gm;
        std::vector< std::size_t > len_pm;

        std::vector< std::map< std::size_t, double >> vec_len_gmpm( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_gm  ( bed_samples.size(), std::map< std::size_t, double >() );
        std::vector< std::map< std::size_t, double >> vec_len_pm  ( bed_samples.size(), std::map< std::size_t, double >() );

        std::vector< std::pair< double, std::string >> temp_gmpm;
        std::vector< std::pair< double, std::string >> temp_gm;
        std::vector< std::pair< double, std::string >> temp_pm;

        std::tuple< double, std::size_t, std::size_t > df_tuple_gmpm;
        std::tuple< double, std::size_t, std::size_t > df_tuple_gm;
        std::tuple< double, std::size_t, std::size_t > df_tuple_pm;

        double pm = 0.0;
        double sum_gmpm = 0.0;
        double sum_gm   = 0.0;
        double sum_pm   = 0.0;

        std::string res_gmpm;
        std::string res_gm;
        std::string res_pm;

        for( auto& anno : ano_len_idx.first )
        {
            vec_gmpm.clear();
            vec_gm  .clear();
            vec_pm  .clear();

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                pm = counting_tables[ smp ].first[ anno ][ 0 ] * counting_tables[ smp ].second[ anno ][ 0 ];

                vec_gmpm.emplace_back( counting_tables[ smp ].first[ anno ][ 0 ] );
                vec_gm  .emplace_back( counting_tables[ smp ].first[ anno ][ 0 ] - pm );
                vec_pm  .emplace_back( pm );

                for( auto& len : ano_len_idx.second )
                {
                    pm = counting_tables[ smp ].first[ anno ][ len ] * counting_tables[ smp ].second[ anno ][ len ];

                    vec_len_gmpm[ smp ][ len ] = counting_tables[ smp ].first[ anno ][ len ];
                    vec_len_gm  [ smp ][ len ] = counting_tables[ smp ].first[ anno ][ len ] - pm;
                    vec_len_pm  [ smp ][ len ] = pm;
                }
            }

            sum_gmpm = get_sum( vec_gmpm );
            sum_gm   = get_sum( vec_gm   );
            sum_pm   = get_sum( vec_pm   );

            df_tuple_gmpm = get_difference( vec_gmpm );
            df_tuple_gm   = get_difference( vec_gm   );
            df_tuple_pm   = get_difference( vec_pm   );

            len_gmpm = get_length_difference( vec_len_gmpm );
            len_gm   = get_length_difference( vec_len_gm   );
            len_pm   = get_length_difference( vec_len_pm   );

            res_gmpm = "\n" + anno + "\t" + std::to_string( sum_gmpm ) + "\t";
            res_gm   = "\n" + anno + "\t" + std::to_string( sum_gm   ) + "\t";
            res_pm   = "\n" + anno + "\t" + std::to_string( sum_pm   ) + "\t";

            if( std::get<0>( df_tuple_gmpm ) > 1 )
                res_gmpm += std::to_string( std::get<0>( df_tuple_gmpm )) + "\t" + bed_samples[ std::get<1>( df_tuple_gmpm )].first + ":" + bed_samples[ std::get<2>( df_tuple_gmpm )].first;
            else
                res_gmpm += "0\t-:-";

            if( std::get<0>( df_tuple_gm ) > 1 )
                res_gm = std::to_string( std::get<0>( df_tuple_gm )) + "\t" + bed_samples[ std::get<1>( df_tuple_gm )].first + ":" + bed_samples[ std::get<2>( df_tuple_gm )].first;
            else
                res_gm = "0\t-:-";

            if( std::get<0>( df_tuple_pm ) > 1 )
                res_pm = std::to_string( std::get<0>( df_tuple_pm )) + "\t" + bed_samples[ std::get<1>( df_tuple_pm )].first + ":" + bed_samples[ std::get<2>( df_tuple_pm )].first;
            else
                res_pm = "0\t-:-";

            res_gmpm += "\t" + std::to_string( get_length_distance( len_gmpm ));
            res_gm   += "\t" + std::to_string( get_length_distance( len_gm   ));
            res_pm   += "\t" + std::to_string( get_length_distance( len_pm   ));

            for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
            {
                res_gmpm += "\t" + std::to_string( len_gmpm[ smp ]);
                res_gm   += "\t" + std::to_string( len_gm  [ smp ]);
                res_pm   += "\t" + std::to_string( len_pm  [ smp ]);
            }

            temp_gmpm.emplace_back( std::make_pair( sum_gmpm, res_gmpm ));
            temp_gm  .emplace_back( std::make_pair( sum_gm  , res_gm   ));
            temp_pm  .emplace_back( std::make_pair( sum_pm  , res_pm   ));
        }

        sort_difference( temp_gmpm );
        sort_difference( temp_gm   );
        sort_difference( temp_pm   );

        for( auto& output : temp_gmpm ) out_gmpm << output.second;
        for( auto& output : temp_gm   ) out_gm   << output.second;
        for( auto& output : temp_pm   ) out_pm   << output.second;

        out_gmpm << "\n";
        out_gm   << "\n";
        out_pm   << "\n";

        out_gmpm.close();
        out_gm.close();
        out_pm.close();
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
