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

    void rerange_un_annotated_start_end(
            std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>& unannos,
            std::map< std::string, std::tuple< std::string, std::size_t, std::size_t >>& un_annotated_map_list,
            const std::size_t& max_anno_merge_size,
            const std::size_t& start,
            const std::size_t& end
            )
    {
        for( auto sit = unannos.find( start ); sit != unannos.end(); ++sit )
            for( auto eit = sit->second.begin(); eit != sit->second.end(); ++eit )
            {
                if( end - start > max_anno_merge_size && end < sit->first ) return;
                for( auto& region_id : eit->second )
                {
                    std::get<1>( un_annotated_map_list[ region_id ] ) = start;
                    std::get<2>( un_annotated_map_list[ region_id ] ) = end;
                }
            }
    }

    void formating_un_annotated_start_end(
            std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>& unannos,
            std::map< std::string, std::tuple< std::string, std::size_t, std::size_t >>& un_annotated_map_list,
            const std::size_t& max_anno_merge_size
            )
    {
        std::size_t start = unannos.begin()->first;
        std::size_t end = unannos.begin()->second.begin()->first;

        for( auto sit = unannos.begin(); sit != unannos.end(); ++sit )
            for( auto eit = sit->second.begin(); eit != sit->second.end(); ++eit )
            {
                if( eit->first - start > max_anno_merge_size && end < sit->first )
                {
                    rerange_un_annotated_start_end(
                            unannos,
                            un_annotated_map_list,
                            max_anno_merge_size,
                            start,
                            end
                            );
                    
                    start = sit->first;
                }

                end = eit->first;
            }

        rerange_un_annotated_start_end(
                unannos,
                un_annotated_map_list,
                max_anno_merge_size,
                start,
                end
                );
    }

    std::vector< std::vector< std::string >> emplace_un_annotated( const auto& anno_rawbed, auto& un_annotated_map_list, auto& un_annotated_checking_map )
    {
        std::vector< std::vector< std::string >> anno_info_temp;
        std::map< std::string, std::vector< std::string >> biotype_map;

        biotype_map[ "un_annotated" ].emplace_back( "un_annotated" );
        biotype_map[ "un_annotated" ].emplace_back(
                anno_rawbed.chromosome_ + ":" +
                std::to_string( anno_rawbed.start_ ) + "-" +
                std::to_string( anno_rawbed.end_ )
                );

        if( un_annotated_map_list.find( biotype_map[ "un_annotated" ][1] ) == un_annotated_map_list.end() )
        {
            un_annotated_map_list[ biotype_map[ "un_annotated" ][1] ] = std::make_tuple( anno_rawbed.chromosome_, anno_rawbed.start_, anno_rawbed.end_ );
            un_annotated_checking_map[ anno_rawbed.chromosome_ ][ anno_rawbed.start_ ][ anno_rawbed.end_ ].emplace_back( biotype_map[ "un_annotated" ][1] );
        }

        anno_info_temp.emplace_back( biotype_map[ "un_annotated" ] );
        return anno_info_temp;
    }

    void filtering(
            std::vector< BedSampleType >& bed_samples,
            std::vector< std::string >& biotype_list,
            const bool& is_keep_other_biotype,
            const std::size_t& max_anno_merge_size
            )
    {
		Filters run_filter;
        std::set< std::string > biotype_set;
        ParaThreadPool smp_parallel_pool( bed_samples.size() );

        for( auto& biotype : biotype_list ) biotype_set.emplace( biotype );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &run_filter, &bed_samples, &biotype_list, &is_keep_other_biotype, &max_anno_merge_size, &biotype_set, this ] () mutable
            {
                std::vector< bool > temp;
                std::vector< std::vector< std::string >> anno_info_temp;
                std::map< std::string, std::vector< std::string >> biotype_map;

                //          regionID                    chr         start       end
                std::map< std::string, std::tuple< std::string, std::size_t, std::size_t >> un_annotated_map_list;
                std::map< std::string, std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>> un_annotated_checking_map;
                //          chr                     start                 end                       regionID

                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    biotype_map.clear();
                    anno_info_temp.clear();

                    anno_rawbed = run_filter.Filter( anno_rawbed );

                    if( anno_rawbed.annotation_info_.empty() || anno_rawbed.annotation_info_[0].empty() )
                    {
                        anno_rawbed.annotation_info_
                            = emplace_un_annotated( anno_rawbed, un_annotated_map_list, un_annotated_checking_map );
                        continue;
                    }

                    for( auto& info : anno_rawbed.annotation_info_ )
                    {
                        for( int i = 0; i < info.size(); i+=2 )
                        {
                            biotype_map[ info[i] ].emplace_back( info[ i   ] );
                            biotype_map[ info[i] ].emplace_back( info[ i+1 ] );
                        }
                    }

                    for( auto& biotpye : biotype_list )
                    {
                        if( biotype_map.find( biotpye ) == biotype_map.end() ) continue;
                        anno_info_temp.emplace_back( biotype_map[ biotpye ] );
                        if( !is_keep_other_biotype ) break;
                    }

                    if( is_keep_other_biotype )
                    {
                        for( auto& biotpye : biotype_map )
                        {
                            if( biotype_set.find( biotpye.first ) != biotype_set.end() ) continue;
                            anno_info_temp.emplace_back( biotpye.second );
                        }
                    }

                    anno_rawbed.annotation_info_ = anno_info_temp;
                }

                for( auto& unannos : un_annotated_checking_map )
                    formating_un_annotated_start_end( unannos.second, un_annotated_map_list, max_anno_merge_size );

                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    if( anno_rawbed.annotation_info_.empty() || anno_rawbed.annotation_info_[0].empty() ) continue;
                    if( anno_rawbed.annotation_info_[0][0] == "un_annotated" )
                    {
                        anno_rawbed.annotation_info_[0][1] =
                            std::get<0>( un_annotated_map_list[ anno_rawbed.annotation_info_[0][1] ]) + ":" +
                            std::to_string( std::get<1>( un_annotated_map_list[ anno_rawbed.annotation_info_[0][1] ])) + "-" +
                            std::to_string( std::get<2>( un_annotated_map_list[ anno_rawbed.annotation_info_[0][1] ])) ;
                    }
                }
            });
        }
        smp_parallel_pool.flush_pool();
    }
};

} // end of namespace algorithm
} // end of namespace ago
