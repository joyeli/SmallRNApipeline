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

    void filtering( std::vector< BedSampleType >& bed_samples, std::vector< std::string >& biotype_list, const bool& is_keep_other_biotype )
    {
		Filters run_filter;
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &run_filter, &bed_samples, &biotype_list, &is_keep_other_biotype ] () mutable
            {
                std::vector< bool > temp;
                std::vector< std::vector< std::string >> anno_info_temp;
                std::map< std::string, std::vector< std::string >> biotype_map;

                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    biotype_map.clear();
                    anno_info_temp.clear();

                    anno_rawbed = run_filter.Filter( anno_rawbed );

                    if( anno_rawbed.annotation_info_.size() == 0 )
                    {
                        biotype_map[ "un_annotated" ].emplace_back( "un_annotated" );
                        biotype_map[ "un_annotated" ].emplace_back(
                            anno_rawbed.chromosome_ + ":" +
                            std::to_string( anno_rawbed.start_ ) + "-" +
                            std::to_string( anno_rawbed.end_ )
                        );

                        anno_info_temp.emplace_back( biotype_map[ "un_annotated" ] );
                        anno_rawbed.annotation_info_ = anno_info_temp;
                        continue;
                    }

                    for( auto& info : anno_rawbed.annotation_info_ )
                    {
                        if( info.size() == 0 )
                        {
                            biotype_map[ "un_annotated" ].emplace_back( "un_annotated" );
                            biotype_map[ "un_annotated" ].emplace_back(
                                anno_rawbed.chromosome_ + ":" +
                                std::to_string( anno_rawbed.start_ ) + "-" +
                                std::to_string( anno_rawbed.end_ )
                            );

                            anno_info_temp.emplace_back( biotype_map[ "un_annotated" ] );
                            anno_rawbed.annotation_info_ = anno_info_temp;
                            continue;
                        }

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

                        // this will keep annotation who is not in the biotype_list
                        if( is_keep_other_biotype )
                            anno_rawbed.annotation_info_ = anno_info_temp;
                    }
                    else anno_info_temp.emplace_back( biotype_map.begin()->second );

                    // this will drop the annotation who is not in the biotype_list
                    if( !is_keep_other_biotype )
                        anno_rawbed.annotation_info_ = anno_info_temp;
                }
            });
        }
        smp_parallel_pool.flush_pool();
    }

};

} // end of namespace algorithm
} // end of namespace ago
