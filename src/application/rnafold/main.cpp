#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

void get_genome( const auto& genome, auto& genome_table )
{
    std::ifstream archive( genome );
    boost::archive::binary_iarchive archive_in( archive );
    archive_in & genome_table;
    archive.close();
}

std::string reverse( std::string seq )
{
    std::transform( seq.begin(), seq.end(), seq.begin(), []( char c )
    {
        switch (c)
        {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'N': c = 'N'; break;
        }

        return c;
    });

    std::reverse( seq.begin(), seq.end() );
    return seq;
}

void get_gtf( const auto& gencode, auto& genome_table, auto& miRs )
{
    std::fstream file( gencode, std::ios::in );

    std::string line;
    std::size_t start;
    std::size_t end;

    std::vector< std::string > split1;
    std::vector< std::string > split2;

    while( std::getline( file, line ))
    {
        if( line.substr( 0, 2 ) == "##" ) continue;
        split1 = explode( line, "\t" );

        if( split1[2] != "gene" ) continue;
        split2 = explode( split1[8], "gene_type \"" );

        if( split2[1].substr( 0, 5 ) != "miRNA" ) continue;
        split2 = explode( split1[8], "gene_name \"" );
        split2 = explode( split2[1], "\"; " );

        start = std::stoi( split1[3] ) -1;
        end   = std::stoi( split1[4] );

        miRs[ ">" + split2[0] + "\t" + std::to_string( start ) + "\t" + split1[4] ]
            = boost::to_upper_copy< std::string >( split1[6] == "-"
            ? reverse( genome_table[ split1[0] ].substr( start, end - start ))
            : genome_table[ split1[0] ].substr( start, end - start ));
    }

    file.close();
}

void get_mirtron( const auto& mirtron, auto& genome_table, auto& miRs  )
{
    std::fstream file( mirtron, std::ios::in );
    std::vector< std::string > split;

    std::string line;
    std::size_t start;
    std::size_t end;

    while( std::getline( file, line ))
    {
        split = explode( line, "\t" );
        start = std::stoi( split[1] );
        end   = std::stoi( split[2] );

        miRs[ ">" + split[4] + "\t" + split[1] + "\t" + split[2] ] = boost::to_upper_copy< std::string >( genome_table[ split[0] ].substr( start, end - start ));
    }

    file.close();
}

std::string exec( const char* cmd )
{
    char buffer[128];
    std::string result = "";

    FILE* pipe = popen( cmd, "r" );
    if( !pipe ) return "ERROR";

    while( !feof( pipe ))
    {
        if( fgets( buffer, 128, pipe ) != NULL )
            result += buffer;
    }

    pclose( pipe );
    return result;
}

int main( int argc, char** argv )
{
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );

    std::string genome, rnafold, mountain, gencode, mirtron;
    std::size_t nthread;

    options.add_options()( "help,h" , "Print help messages" )
        ( "genome,g"  , boost::program_options::value< std::string >( &genome   )->required(), "Set genome archive file" )
        ( "rnafold,r" , boost::program_options::value< std::string >( &rnafold  )->required(), "Set path of RNAfold" )
        ( "mountain,m", boost::program_options::value< std::string >( &mountain )->required(), "Set path of mountain.pl" )
        ( "gencode,c" , boost::program_options::value< std::string >( &gencode  )->required(), "Set input main gtf annotation file from gencode" )
        ( "mirtron,t" , boost::program_options::value< std::string >( &mirtron  )->default_value(""),"Set input mirtron bed annotation file from liftover bed format" )
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

    std::map< std::string, std::string > genome_table;
    std::map< std::string, std::string > miRs;
    std::map< std::string, std::string > miRs_folds;
    std::map< std::string, std::string > miRs_mfes;

    std::map< std::string, std::vector< std::string >> miRs_enps;

    std::vector< std::string > split_res1;
    std::vector< std::string > split_res2;
    std::vector< std::string > split_colu;

    std::string results;
    std::string this_mirna;

    std::size_t line_count;
    std::size_t and_count;

    get_genome( genome, genome_table );
    get_gtf( gencode, genome_table, miRs );
    get_mirtron( mirtron, genome_table, miRs );

    std::ofstream output( "./miRs.fa" );

    for( auto& miR : miRs )
        output << miR.first << "\n" << miR.second << "\n";

    output.close();

    results = exec(( rnafold + " -p --noLP < ./miRs.fa" ).c_str());
    split_res1 = explode( results, "\n" );

    for( auto& line : split_res1 )
    {
        if( line == "" ) continue;

        if( line_count == 2 )
        {
            split_colu = explode( line, " " );

            miRs_folds[ this_mirna ] = split_colu[0];
            miRs_mfes[  this_mirna ] = split_colu[1].substr( 1, split_colu[1].length() -2 );
        }
        
        if( line.at(0) == '>' )
        {
            line_count = 0;
            this_mirna = line;

            split_colu = explode( line, "\t" );
            results = exec(( mountain + " ./" + split_colu[0].substr( 1, split_colu[0].length() -1 ) + "_dp.ps" ).c_str()); 

            split_res2 = explode( results, "\n" );
            and_count = 0;

            for( auto& value : split_res2 )
            {
                if( value == "" ) continue;
                if( value == "&" ) and_count++;
                if( and_count == 2 )
                {
                    if( value == "&" ) continue;
                    split_colu = explode( value, " " );
                    miRs_enps[ this_mirna ].emplace_back( split_colu[ split_colu.size() -1 ]);
                }
            }
        }

        line_count++;
    }

    output.open( "./miRNA_RNAFold.tsv" );

    for( auto& miR : miRs )
    {
        if( miRs_enps[ miR.first ].empty() ) continue;

        output
            << miR.first.substr( 1, miR.first.length() -1 ) << "\t"
            << miRs_mfes[ miR.first ] << "\t" << miR.second << "\t"
            << miRs_folds[ miR.first ];

        for( auto& value : miRs_enps[ miR.first ] )
            output << "\t" << value;

        output << "\n";
    }

    output.close();
    return 0;
}
