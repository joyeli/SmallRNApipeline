#pragma once
#include <AGO/format/md_sam.hpp>
#include <pokemon/format/raw_bed.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

namespace ago {
namespace format {

// the start from RawBedBase is 1-base
struct MDRawBed : public RawBedBase
{
    char strand_;
    uint64_t end_;
    size_t is_filtered_;

    double ppm_;

    std::string chromosome_;
    std::vector< std::vector< std::string >> annotation_info_;
    std::vector< std::vector< bool >> is_on_biotype_list_;

    std::map< std::size_t, char > md_map; // map< index, nt > for mismatch
    std::set< std::size_t > tc_set; // set< index > fot T-to-C
    std::string seed_md_tag;
    std::string tail;

    MDRawBed( void )
        : RawBedBase()
    {}

    MDRawBed( const MDSam<>& samin )
        : RawBedBase( samin )
        , strand_( get_strand( samin ))
        , end_( get_end( samin ))
        , chromosome_( get_chromosome( samin ))
        , is_filtered_(0)
        , ppm_(0)
        , annotation_info_(0)
        , is_on_biotype_list_(0)
        , md_map( samin.md_map )
        , tc_set( samin.tc_set )
        , seed_md_tag( get_md_in_seed( samin ))
        , tail( std::get<2>( std::get<11>( samin.data )))
    {
        this->tail_length_ = (uint8_t)( tail.length() );
    }

    char get_strand( const MDSam<>& samin )
    {
        if( std::get<1>( samin.data ) == 0 )
        {
            return '+';
        }
        else
        {
            return '-';
        }
    }

    std::size_t get_genome_match_length( const MDSam<>& samin )
    {
        std::size_t gm_length = 0;
        std::vector< std::string > cigar;
    
        if( std::get<1>( samin.data ) == 16 )
        {
            boost::split( cigar, std::get<5>( samin.data ), boost::is_any_of( "S" ));
            switch( cigar.size() )
            {
                case 1 : gm_length = std::stoi( cigar[0].substr( 0, cigar[0].length() -1 )); break;
                default: gm_length = std::stoi( cigar[1].substr( 0, cigar[1].length() -1 )); break;
            }
        }
        else
        {
            boost::split( cigar, std::get<5>( samin.data ), boost::is_any_of( "M" ));
            gm_length = std::stoi( cigar[0] );
        }
        return gm_length;
    }

    uint64_t get_end( const MDSam<>& samin )
    {
        return std::get<3>( samin.data ) + get_genome_match_length( samin );
    }

    std::string get_chromosome( const MDSam<>& samin )
    {
        if( std::get<2>( samin.data ).size() > 3 && std::get<2>( samin.data ).substr( 0, 3 ) == "chr" )
        {
            return std::get<2>( samin.data );
        }
        else
        {
            return "chr" + std::get<2>( samin.data );
        }
    }

    std::string get_genome_match_raw_seed_sequence( const MDSam<>& samin )
    {
        std::string seq = std::get<9>( samin.data );
        if( std::get<1>( samin.data ) == 16 )
        {
            std::reverse( seq.begin(), seq.end() );
			std::transform( seq.begin(), seq.end(), seq.begin(),
                    [ this ]( char c ){ return this->complement( c ); });
        }
        return seq.substr( 1, 7 );
    }

    std::string get_md_in_seed( const MDSam<>& samin )
    {
        std::string tag = "";
        std::string seed = get_genome_match_raw_seed_sequence( samin );
        for( auto& md : md_map )
        {
            if( md.first > 0 && md.first < 8 )
            {
                tag += std::to_string( md.first -1 ) + seed.at( md.first -1 );
            }
        }
        return tag;
    }

	std::string getTail()
	{
		return this->tail;
	}

	char complement( char c )
	{
		switch( c )
        {
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
		}
		return c;
	}

	std::string getReadSeq( std::map< std::string, std::string >& gGenometable )
	{
		std::string read_seq = gGenometable[ chromosome_ ].substr( this->start_ -1, this->length_ - this->tail_length_ );
		std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper );
		if( strand_ == '-' )
		{
			std::reverse( read_seq.begin(), read_seq.end() );
			std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(),
                    [ this ]( char c ){ return this->complement( c ); });
		}
        read_seq += getTail();
		return std::move( read_seq );
	}

    std::string getMD()
    {
        std::string md = "";
        for( auto& m : md_map )
            md += std::to_string( m.first ) + m.second;
        return md;
    }

    std::string getTC()
    {
        std::string tc = "";
        for( auto& t : tc_set )
            tc += std::to_string( t ) + "C";
        return tc;
    }

	friend class boost::serialization::access;
	template< class Archive >
	void serialize( Archive &ar, const unsigned int version )
	{  
        ar & this->chr_idx_; 
		ar & this->tail_length_; 
		ar & this->length_; 
		ar & this->multiple_alignment_site_count_; 
		ar & this->reads_count_; 
		ar & this->start_; 
		ar & this->tail_mismatch_;

		ar & strand_;
		ar & end_;
		ar & is_filtered_;
		ar & ppm_;
		ar & chromosome_;
		ar & annotation_info_;
		ar & is_on_biotype_list_;

        ar & md_map;
        ar & tc_set;
        ar & seed_md_tag;
        ar & tail;
	}
};

} // format
} // ago
