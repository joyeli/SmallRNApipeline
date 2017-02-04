#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

void emplace_back( std::vector< std::vector< std::string >>& bed_vec, const std::vector< std::string >& split_vec, const std::string& type )
{
    uint64_t start = (uint64_t)(std::stol( split_vec[3] )) -1;

    std::vector< std::string > split_type_vec = explode( split_vec[8], type + "_type \"" );
    std::vector< std::string > split_type     = explode( split_type_vec[1], "\"; " );
    std::vector< std::string > split_name_vec = explode( split_vec[8], type + "_name \"" );
    std::vector< std::string > split_name     = explode( split_name_vec[1], "\"; " );

    bed_vec.emplace_back(
        std::vector< std::string >{
              split_vec[0]
            , std::to_string( start )
            , split_vec[4]
            , split_vec[6]
            , split_type[0]
            , split_name[0]
        }
    );
}

void emplace_bed( std::vector< std::vector< std::string >>& res_vec, const std::vector< std::vector< std::string >>& bed_vec )
{
    res_vec.emplace_back( bed_vec[0] );

    for( size_t i = 1; i < bed_vec.size(); ++i )
    {
        if( bed_vec[0][4] != bed_vec[i][4] )
        {
            res_vec.emplace_back( bed_vec[i] );
        }
    }
}

bool is_filter_type( const std::string& type )
{
    if(
           type == "vaultRNA"
        || type == "nonsense_mediated_decay"
        || type == "non_stop_decay"
        || type == "retained_intron"
    )
    {
        return true;
    }

    return false;
}

int main( int argc, char** argv )
{
    switch( argc )
    {
        case 2  : break;
        default : throw std::runtime_error( "This version of converter is base on gencode v25/vM12\n  Usage:  ./EXE gencode.gtf" );
    }

    std::string line;
    std::string name = argv[1];

    std::vector< std::string > split_vec;
    std::vector< std::vector< std::string >> bed_vec;
    std::vector< std::vector< std::string >> res_vec;

    std::fstream file( name, std::ios::in );

    while( std::getline( file, line ))
    {
        if( line.substr(0,1) == "#" )
            continue;

        split_vec = explode( line, "\t" );

        if( split_vec[2] == "transcript" )
        {
            emplace_back( bed_vec, split_vec, "transcript" );
            continue;
        }

        if( split_vec[2] == "gene" )
        {
            if( !bed_vec.empty() )
            {
                emplace_bed( res_vec, bed_vec );
            }

            bed_vec.clear();
            emplace_back( bed_vec, split_vec, "gene" );
        }
    }

    file.close();
    emplace_bed( res_vec, bed_vec );

    std::ofstream output( name + ".bed" );

    for( auto& bed : res_vec )
    {
        if( !is_filter_type( bed[4] ))
        {
            output << bed[0] << "\t" << bed[1] << "\t" << bed[2] << "\t" << bed[3] << "\t" << bed[4] << "\t" << bed[5] << "\n";
        }
    }

    output.close();

    return 0;
}

