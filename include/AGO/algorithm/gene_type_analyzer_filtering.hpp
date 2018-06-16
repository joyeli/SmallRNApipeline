#pragma once
#include <AGO/format/md_rawbed.hpp>
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerFiltering
{
  public:

    GeneTypeAnalyzerFiltering()
    {}

    using DropTypeList = boost::mpl::vector<
        boost::mpl::vector< boost::mpl::string< 'rm', 'sk' >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>>;

    // using Filters = FilterWorker< AnnotationRawBed<>, DropTypeList >;
    using Filters = FilterWorker< ago::format::MDRawBed, DropTypeList >;

    void filtering( std::vector< BedSampleType >& bed_samples, std::vector< std::string >& biotype_list, bool& is_filter_drop )
    {
		Filters run_filter;
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        std::set< std::string > biotype_set;
        // for( auto& biotype : biotype_list ) biotype_set.emplace( biotype );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &run_filter, &bed_samples, &biotype_set, &is_filter_drop, &biotype_list ] () mutable
            {
                std::vector< bool > temp;
                std::vector< std::vector< std::string >> anno_info_temp;
                std::map< std::string, std::vector< std::string >> biotype_map;

                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    anno_rawbed = run_filter.Filter( anno_rawbed );

                    for( auto& info : anno_rawbed.annotation_info_ )
                    {
                        for( int i = 0; i < info.size(); i+=2 )
                        {
                            biotype_map[ info[i] ].emplace_back( info[ i   ] );
                            biotype_map[ info[i] ].emplace_back( info[ i+1 ] );
                        }
                    }

                    if( biotype_map.size() > 1 )
                    {
                        for( auto& biotpye : biotype_list )
                        {
                            if( biotype_map.find( biotpye ) != biotype_map.end() )
                            {
                                anno_info_temp.emplace_back( biotype_map[ biotpye ] );
                            }
                        }

                        anno_rawbed.annotation_info_ = anno_info_temp;
                    }

                    anno_info_temp.clear();
                    biotype_map.clear();

                    // this is for anno_rawbed.is_on_biotype_list_
                    // for( auto& info : anno_rawbed.annotation_info_ )
                    // {
                    //     for( int i = 0; i < info.size(); i+=2 )
                    //     {
                    //         if( biotype_set.find( info[i] ) != biotype_set.end() )
                    //         {
                    //             temp.emplace_back( true );
                    //             temp.emplace_back( true );
                    //         }
                    //         else
                    //         {
                    //             temp.emplace_back( false );
                    //             temp.emplace_back( false );
                    //         }
                    //     }

                    //     anno_rawbed.is_on_biotype_list_.emplace_back( temp );
                    //     temp.clear();
                    // }
                }
            });
        }
        smp_parallel_pool.flush_pool();
    }

};

} // end of namespace algorithm
} // end of namespace ago
