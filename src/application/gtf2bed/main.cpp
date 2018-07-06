#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

std::vector< std::vector< std::string >> read_gtf( const std::string& anno, const bool& main_annotaion = false )
{
    std::fstream file( anno, std::ios::in );
    std::string line;

    std::vector< std::vector< std::string >> bed_vec;
    std::vector< std::string > split_vec;

    std::vector< std::string > split_type_vec;
    std::vector< std::string > split_type;

    std::vector< std::string > split_name_vec;
    std::vector< std::string > split_name;

    while( std::getline( file, line ))
    {
        if( line.substr( 0, 2 ) == "##" )
            continue;

        split_vec = explode( line, "\t" );
        
        if( main_annotaion == true && split_vec[2] != "gene" )
            continue;

        split_type_vec = explode( split_vec[8], "gene_type \"" );
        split_type     = explode( split_type_vec[1], "\"; " );

        split_name_vec = explode( split_vec[8], "gene_name \"" );
        split_name     = explode( split_name_vec[1], "\"; " );

        bed_vec.emplace_back( std::vector< std::string >{
                  split_vec[0]
                , std::to_string( std::stol( split_vec[3] ) -1 )
                , split_vec[4]
                , split_vec[6]
                , split_type[0]
                , split_name[0]
                });
    }

    file.close();
    return bed_vec;
}

std::vector< std::vector< std::string >> read_bed6( const std::string& rmsk, const std::string& type )
{
    std::fstream file( rmsk, std::ios::in );
    std::string line;

    std::vector< std::vector< std::string >> bed_vec;
    std::vector< std::string > split_vec;

    while( std::getline( file, line ))
    {
        split_vec = explode( line, "\t" );
        bed_vec.emplace_back( std::vector< std::string >{
                  split_vec[0]
                , split_vec[1]
                , split_vec[2]
                , split_vec[5]
                , type
                , split_vec[3] + "_" + split_vec[4]
                });
    }

    file.close();
    return bed_vec;
}

std::vector< std::vector< std::string >> read_mirtron( const std::string& mtro )
{
    std::fstream file( mtro, std::ios::in );
    std::string line;

    std::vector< std::vector< std::string >> bed_vec;
    std::vector< std::string > split_vec;

    while( std::getline( file, line ))
    {
        split_vec = explode( line, "\t" );
        bed_vec.emplace_back( std::vector< std::string >{
                  split_vec[0]
                , split_vec[1]
                , split_vec[2]
                , split_vec[3]
                , "mirtron"
                , split_vec[4]
                });
    }

    file.close();
    return bed_vec;
}

std::map< std::string, std::string > get_options( int& argc, char** argv )
{
    std::map< std::string, std::string > args;
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );
    options.add_options()( "help,h" , "Print help messages" )
        ( "grch,g", boost::program_options::value< std::string >()->required(),"Set human genome version" )
        ( "name,n", boost::program_options::value< std::string >()->required(),"Set output file name" )
        ( "anno,a", boost::program_options::value< std::string >()->required(),"Set input main gtf annotation file from Gencode" )
        ( "rmsk,r", boost::program_options::value< std::string >()->default_value(""),"Set input rmsk bed6 annotation file from UCSC" )
        ( "mtro,m", boost::program_options::value< std::string >()->default_value(""),"Set input mirtron bed annotation file from liftover bed format" )
        ( "trna,t", boost::program_options::value< std::string >()->default_value(""),"Set input tRNA gtf annotation file from GenCode" )
        ( "plya,p", boost::program_options::value< std::string >()->default_value(""),"Set input polyA gtf annotation file from GenCode" )
        ( "sudo,s", boost::program_options::value< std::string >()->default_value(""),"Set input pseudogene gtf annotation file from GenCode" )
        ;
    try {
        boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), op );
        if( op.count( "help" )) {
            std::cout << "\n" << options << "\n";
            exit(0);
        }
        boost::program_options::notify( op );
    }
    catch( boost::program_options::error& error ) {
        std::cerr << "\nERROR: " << error.what() << "\n" << options << "\n";
        exit(1);
    }
    for( auto& arg : op ) {
        args[ arg.first ] = arg.second.as< std::string >();
    }
    return args;
}

int main( int argc, char** argv )
{
    std::map< std::string, std::string > args( get_options( argc, argv ));
    std::map< std::string, std::size_t > counting_type;
    
    auto& grch = args[ "grch" ];
    auto& name = args[ "name" ];
    auto& anno = args[ "anno" ];
    auto& rmsk = args[ "rmsk" ];
    auto& mtro = args[ "mtro" ];
    auto& trna = args[ "trna" ];
    auto& plya = args[ "plya" ];
    auto& sudo = args[ "sudo" ];

    std::vector< std::string > trna_check;
    std::vector< std::vector< std::string >> bed_vec;
    std::vector< std::vector< std::string >> anno_bed_vec;
    std::vector< std::vector< std::string >> mirtron_bed_vec;

    bool is_mirtron = false;

    anno_bed_vec = read_gtf( anno, true );
    mirtron_bed_vec = read_mirtron( mtro );

    for( auto& bed : anno_bed_vec )
    {
        if( bed[4] == "miRNA" )
        {
            for( auto& mirtron : mirtron_bed_vec )
            {
                if( bed[0] == mirtron[0] && bed[3] == mirtron[3]   &&
                    std::stoi( bed[2] ) >= std::stoi( mirtron[1] ) &&
                    std::stoi( bed[1] ) <= std::stoi( mirtron[2] ))
                {
                    is_mirtron = true;
                    break;
                }
            }

            if( is_mirtron )
            {
                is_mirtron = false;
                continue;
            }

            // bed_vec.emplace_back( bed );
            // bed[4] = "miRNA_mirtron";
        }
        
        bed_vec.emplace_back( bed );
    }

    for( auto& bed : mirtron_bed_vec )
    {
        // bed_vec.emplace_back( bed );
        // bed[4] = "miRNA_mirtron";
        bed_vec.emplace_back( bed );
    }

    anno_bed_vec.clear();

    anno_bed_vec = read_bed6( rmsk, "rmsk" );
    bed_vec.insert( bed_vec.end(), anno_bed_vec.begin(), anno_bed_vec.end() );
    anno_bed_vec.clear();

    if( trna != "" )
    {
        anno_bed_vec = read_gtf( trna );
        bed_vec.insert( bed_vec.end(), anno_bed_vec.begin(), anno_bed_vec.end() );
        anno_bed_vec.clear();
    }

    if( plya != "" )
    {
        anno_bed_vec = read_gtf( plya );
        bed_vec.insert( bed_vec.end(), anno_bed_vec.begin(), anno_bed_vec.end() );
        anno_bed_vec.clear();
    }

    if( sudo != "" )
    {
        anno_bed_vec = read_gtf( sudo );
        bed_vec.insert( bed_vec.end(), anno_bed_vec.begin(), anno_bed_vec.end() );
        anno_bed_vec.clear();
    }

    std::ofstream output;
    output.open( name + "." + grch + ".bed" );

    for( auto& bed : bed_vec )
    {
        trna_check = explode( bed[4], "tRNA" );

        if( trna_check.size() != 1 )
        {
            bed[5] = bed[4] + "_" + bed[5];
            bed[4] = "tRNA";
        }

        output
            << bed[0] << "\t"   // chr
            << bed[1] << "\t"   // start
            << bed[2] << "\t"   // end
            << bed[3] << "\t"   // strand
            << bed[4] << "\t"   // annotation type
            << bed[5] << "\n";  // annotation name

        if( counting_type.find( bed[4] ) != counting_type.end() )
             counting_type[ bed[4] ] += 1;
        else counting_type[ bed[4] ]  = 1;
    }

    output.close();
    output.open( name + "." + grch + ".type" );

    for( auto& type : counting_type )
    {
        output << type.first << "\t" << type.second << "\n";
    }

    output.close();
    return 0;
}
