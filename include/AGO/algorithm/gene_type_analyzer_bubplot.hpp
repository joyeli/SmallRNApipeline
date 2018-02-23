#pragma once
#include <unistd.h>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBubplot
{
    using ChrRangeType = std::tuple< std::string, std::size_t, std::size_t, char >;

  public:

    GeneTypeAnalyzerBubplot()
    {}

    static void output_bubplot(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number,
            auto& genome_table
            )
    {
        for( auto& smp : bed_samples ) if( !boost::filesystem::exists( output_name + smp.first + ".tsv" ))
            boost::filesystem::create_symlink(( "../LenDist/" + smp.first + ".tsv" ).c_str(), ( output_name + smp.first + ".tsv" ).c_str() );

        output_anno_sequence( output_name, bed_samples, biotype, thread_number, genome_table );
    }

    static void output_anno_sequence(
            const std::string& output_name,
            const std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number,
            auto& genome_table
            )
    {
        std::map< std::string, ChrRangeType > chr_mapping = get_chrmap_table( bed_samples, biotype, thread_number );
        std::vector< std::string > out_vec = sequence_formating( chr_mapping, genome_table );

        std::ofstream output( output_name + "anno_seq.tsv" );
        output << "Annotation\t5P\t3P";

        for( auto& res : out_vec )
            output << res;

        output << "\n";
        output.close();
    }

    static std::map< std::string, ChrRangeType > get_chrmap_table(
            const std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            const std::size_t& thread_number
            )
    {
        std::map< std::string, std::vector< std::pair< ChrRangeType, std::size_t >>> chr_mappings;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                for( auto& raw_bed_info : raw_bed.annotation_info_ )
                {
                    for( std::size_t i = 0; i < raw_bed_info.size(); i+=2 )
                    {
                        if( raw_bed_info[i] != biotype ) continue;

                        if( chr_mappings.find( raw_bed_info[ i+1 ] ) == chr_mappings.end() )
                            chr_mappings[ raw_bed_info[ i+1 ]] = std::vector< std::pair< ChrRangeType, std::size_t >>();

                        chr_mappings[ raw_bed_info[ i+1 ]].emplace_back(
                                std::make_pair( std::make_tuple(
                                    raw_bed.chromosome_,
                                    raw_bed.start_,
                                    raw_bed.end_,
                                    raw_bed.strand_ ),
                                    0
                                ));
                    }
                }
            }
        }

        ParaThreadPool parallel_pool( thread_number );
        std::vector< std::pair< std::string, std::vector< std::pair< ChrRangeType, std::size_t >>>> parallel_vec;
        std::vector< std::size_t > parallel_indx;

        std::size_t task_number = thread_number * 10;

        for( auto& anno : chr_mappings )
            parallel_vec.emplace_back( anno );

        chr_mappings.clear();

        for( std::size_t anno = 0; anno < parallel_vec.size(); ++anno )
        {
            parallel_indx.emplace_back( anno );

            if( parallel_indx.size() >= task_number )
                parallel_pool.job_post([ parallel_indx, &parallel_vec ] ()
                {
                    for( auto& idx : parallel_indx )
                    {
                        recursive_merge( parallel_vec[ idx ].second, 0 );
                        range_counting_sort( parallel_vec[ idx ].second );
                    }
                });

            parallel_indx.clear();
        }

        parallel_pool.flush_pool();

        std::map< std::string, ChrRangeType > chr_mapping_res;

        for( auto& anno : parallel_vec )
            chr_mapping_res[ anno.first ] = anno.second[0].first;

        return chr_mapping_res;
    }

    static void recursive_merge( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges, std::size_t start_idx )
    {
        if( start_idx == ranges.size() ) return;
        std::size_t counts  = ranges[ start_idx ].second;

        for( std::size_t idx = 0; idx < ranges.size(); ++idx )
        {
            if( idx == start_idx ) continue;
            
            if( std::get<0>( ranges[ idx ].first ) != std::get<0>( ranges[ start_idx ].first )) continue;
            if( std::get<3>( ranges[ idx ].first ) != std::get<3>( ranges[ start_idx ].first )) continue;

            if( std::get<1>( ranges[ idx ].first ) <= std::get<2>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) >= std::get<1>( ranges[ start_idx ].first ) || 
                std::get<1>( ranges[ idx ].first ) >= std::get<1>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) <= std::get<2>( ranges[ start_idx ].first ) || 
                std::get<1>( ranges[ idx ].first ) <= std::get<1>( ranges[ start_idx ].first ) && std::get<2>( ranges[ idx ].first ) >= std::get<2>( ranges[ start_idx ].first ) )
            {
                if( std::get<1>( ranges[ idx ].first ) < std::get<1>( ranges[ start_idx ].first )) std::get<1>( ranges[ start_idx ].first ) = std::get<1>( ranges[ idx ].first );
                if( std::get<2>( ranges[ idx ].first ) > std::get<2>( ranges[ start_idx ].first )) std::get<2>( ranges[ start_idx ].first ) = std::get<2>( ranges[ idx ].first );
                if( idx < start_idx ) start_idx--;

                ranges.erase( ranges.begin(), ranges.begin() + idx + 1 );
                counts++;
                idx--;
            }
        }

        if( counts == ranges[ start_idx ].second ) recursive_merge( ranges, start_idx +1 );
    }

    static void range_counting_sort( std::vector< std::pair< ChrRangeType, std::size_t >>& ranges )
    {
        std::sort( ranges.begin(), ranges.end(),
            []( const std::pair< ChrRangeType, std::size_t >& a, const std::pair< ChrRangeType, std::size_t >& b )
            { return a.second > b.second; });
    }

    static std::vector< std::string > sequence_formating(
            const std::map< std::string, ChrRangeType >& chr_mapping,
            auto& genome_table
            )
    {
        std::vector< std::string > res_vec;
        std::map< std::string, std::tuple< std::string, std::string >> temp_map;

        for( auto& anno : chr_mapping )
        {
            if( temp_map.find( anno.first.substr( 0, anno.first.length() -3 )) == temp_map.end() )
                temp_map[ anno.first.substr( 0, anno.first.length() -3 )] = std::tuple< std::string, std::string >();

            switch( anno.first.at( anno.first.length() -2 ))
            {
                case '5' : std::get<0>( temp_map[ anno.first.substr( 0, anno.first.length() -3 )]) = get_sequence( anno.second, genome_table ); break;
                case '3' : std::get<1>( temp_map[ anno.first.substr( 0, anno.first.length() -3 )]) = get_sequence( anno.second, genome_table ); break;
            }
        }

        for( auto& anno : temp_map )
        {
            res_vec.emplace_back(
                    "\n" + anno.first +
                    "\t" + std::get<0>( anno.second ) +
                    "\t" + std::get<1>( anno.second ) );
        }

        return res_vec;
    }

    static std::string get_sequence( const ChrRangeType& range, auto& genome_table )
    {
		std::string read_seq =
            genome_table[ std::get<0>( range )].substr( std::get<1>( range ) - 1, std::get<2>( range ) - std::get<1>( range ));

		std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper );

		if( std::get<3>( range ) == '-' )
		{
			std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), []( char c ){ return complement( c ); });
			std::reverse( read_seq.begin(), read_seq.end() );
		}

        return read_seq;
    }

	static char complement( char c )
	{
		switch (c) {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
		}
		return c;
	}
};

} // end of namespace algorithm
} // end of namespace ago
