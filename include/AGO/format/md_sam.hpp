#pragma once
#include <pokemon/format/sam.hpp>

namespace ago {
namespace format {

typedef UserDefineTags< // TAG, 11
    boost::mpl::string<'NH'>, 
    boost::mpl::string<'MD'>, 
    boost::mpl::string<'TL'>
> NHMDTLTags;

typedef std::tuple <
        std::string,                // QNAME    , 0
        int,                        // FLAG     , 1
        std::string,                // RNAME    , 2
        uint64_t,                   // POS      , 3
        int,                        // MAPQ     , 4
        std::string,                // CIGAR    , 5
        std::string,                // RNEXT    , 6
        uint64_t,                   // PNEXT    , 7
        int64_t,                    // TLEN     , 8
        std::string,                // SEQ      , 9
        std::string,                // QUAL     , 10
        NHMDTLTags                  // TAG      , 11
> MDSamTupleType;

template< class TupleType = MDSamTupleType >
struct MDSam : public Sam< TupleType >
{
    std::map< std::size_t, char > md_map; // map< index, nt > for mismatch
    std::set< std::size_t > tc_set; // set< index > fot T-to-C

    MDSam( void )
        : Sam< TupleType >()
    {}

    MDSam( const MDSam< TupleType >& md_sam )
        : Sam< TupleType >( md_sam )
        , md_map( md_sam.md_map )
        , tc_set( md_sam.tc_set )
    {}

	MDSam< TupleType >& operator= ( const MDSam< TupleType >& other )
	{
		this->eof_flag = other.eof_flag;
		this->data = other.data;
        md_map = other.md_map;
        tc_set = other.tc_set;
		return *this;
	}

	MDSam( TupleType& data_in )
		: Sam< TupleType >( data_in ) 
	{
        get_md_tc( this->data );
	}

	MDSam( TupleType&& data_in) 
		: Sam< TupleType >( std::move( data_in ))
	{
        get_md_tc( this->data );
	}

	MDSam( MDSam< TupleType >&& other ) 
        : Sam< TupleType >( std::move( other.data ))
        , md_map( std::move( other.md_map ))
        , tc_set( std::move( other.tc_set ))
	{}

	MDSam< TupleType >& operator= ( MDSam< TupleType > && other )
	{
		this->data = std::move( other.data );
		this->eof_flag = std::move( other.eof_flag );
        md_map = std::move( other.md_map );
        tc_set = std::move( other.tc_set );
		return *this;
	}

	MDSam( bool EofFlag ) 
		: Sam< TupleType >( EofFlag ) 
	{}

	MDSam( std::string& line )
        : Sam< TupleType >( line )
	{
        get_md_tc( this->data, false );
	}

	MDSam( std::string&& line )
        : Sam< TupleType >( line )
	{
        get_md_tc( this->data, false );
	}

	MDSam( std::string& line, bool& allow_t2c )
        : Sam< TupleType >( line )
	{
        get_md_tc( this->data, allow_t2c );
	}

	MDSam( std::string&& line, bool& allow_t2c )
        : Sam< TupleType >( line )
	{
        get_md_tc( this->data, allow_t2c );
	}

    void get_md_tc( TupleType& tdata, bool& allow_t2c )
    {
        auto& flag = std::get<1>( tdata );
        auto& md_tag = std::get<1>( std::get<11>( tdata ));

        std::size_t gm = get_genome_match_length( tdata );
        std::string seq = get_genome_match_sequence( tdata, gm );
        std::string md_str = "";

        std::size_t md_pos = 0;
        std::size_t last_md = 0;

        int last_md_pos = -1;
        bool is_atcg = false;
        bool is_t2c = false;

        for( std::size_t i = 0; i < md_tag.length(); ++i )
        {
            switch( md_tag.at(i) )
            {
                case 'A' : is_atcg = true; break;
                case 'T' : is_atcg = true; break;
                case 'C' : is_atcg = true; break;
                case 'G' : is_atcg = true; break;
                default  : is_atcg = false;
            }
            if( is_atcg )
            {
                for( std::size_t j = last_md_pos +1; j < i; ++j )
                    md_str += md_tag.at(j);

                md_pos = last_md + std::stoi( md_str );

                if( flag == 0  && md_tag.at(i) == 'T' && seq.at( md_pos ) == 'C' ) is_t2c = true;
                if( flag == 16 && md_tag.at(i) == 'A' && seq.at( md_pos ) == 'G' ) is_t2c = true;

                if( !is_t2c || !allow_t2c )
                {
                    if( flag == 0 ) md_map[ md_pos ] = seq.at( md_pos );
                    else switch( seq.at( md_pos ))
                    {
                        case 'A' : md_map[ gm -1 - md_pos ] = 'T'; break;
                        case 'T' : md_map[ gm -1 - md_pos ] = 'A'; break;
                        case 'C' : md_map[ gm -1 - md_pos ] = 'G'; break;
                        case 'G' : md_map[ gm -1 - md_pos ] = 'C'; break;
                    }
                }
                else switch( flag )
                {
                    case 0 : tc_set.emplace( md_pos ); break;
                    default: tc_set.emplace( gm -1 - md_pos );
                }

                is_t2c = false;
                last_md_pos = i;
                last_md = md_pos +1;
                md_str = "";
            }
        }
    }

    std::size_t get_genome_match_length( const TupleType& tdata )
    {
        std::size_t gm_length = 0;
        std::vector< std::string > cigar;
    
        if( std::get<1>( tdata ) == 16 )
        {
            boost::split( cigar, std::get<5>( tdata ), boost::is_any_of( "S" ));
            switch( cigar.size() )
            {
                case 1 : gm_length = std::stoi( cigar[0].substr( 0, cigar[0].length() -1 )); break;
                default: gm_length = std::stoi( cigar[1].substr( 0, cigar[1].length() -1 )); break;
            }
        }
        else
        {
            boost::split( cigar, std::get<5>( tdata ), boost::is_any_of( "M" ));
            gm_length = std::stoi( cigar[0] );
        }
        return gm_length;
    }

    std::string get_genome_match_sequence( const TupleType& tdata, const std::size_t& gm )
    {
        std::string seq = std::get<9>( tdata );
        if( std::get<1>( tdata ) == 16 )
        {
            std::reverse( seq.begin(), seq.end() );
            seq = seq.substr( 0, gm );
            std::reverse( seq.begin(), seq.end() );
        }
        else seq = seq.substr( 0, gm );
        return seq;
    }

	friend class boost::serialization::access;
	template< class Archive >
	void serialize( Archive &ar, const unsigned int version )
	{
		TupleUtility < TupleType, std::tuple_size< TupleType >::value > :: SerializeTuple( this->data, ar, version );
		ar & this->eof_flag;
        ar & md_map;
        ar & tc_set;
	}
};

} // format
} // ago
