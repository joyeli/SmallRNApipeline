#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <Tailor/tailer.hpp>
#include <boost/program_options.hpp>

std::map< std::string, std::string > get_options( int& argc, char** argv );

int main( int argc, char** argv )
{
    std::map< std::string, std::string > args( get_options( argc, argv ));

    std::string readsPath = args[ "fastq" ];
    std::string idxPrefix = args[ "index" ];
    std::string genome    = args[ "genome" ];

    int  threadNum = std::stoi( args[ "thread" ] );
    int  minLen    = std::stoi( args[ "minLen" ] );
    bool allowMM   = std::stoi( args[ "mismatch" ]);

    std::stringstream ssMsg;

    if( genome != "" )
        tailor::buildBWT2( genome, idxPrefix );

    ABWT_table abwtt = tailor::loadBWT2( idxPrefix, &ssMsg );
    // tailor::searchBWT_tail2( std::move( abwtt ), readsPath, threadNum, &ssMsg, minLen, allowMM );

    boost::thread* threads [ threadNum ];
	std::ifstream in{ readsPath };

    std::string s;
    std::stringstream ss;
    bool new_line_check = false;

    while( std::getline( in, s ))
    {
        if( new_line_check )
        {
            ss << "\n";
        }
        else
        {
            new_line_check = true;
        }
        ss << s;
    }

    std::cerr << "\n=================\n" << ss.str() << "\n=================\n";

    for( size_t i = 0 ; i < threadNum; ++i )
    {
        threads[i] = new boost::thread
        {
            tailor::ABWT_threads< ABWT_table >{
                abwtt, &ss, &ssMsg, minLen, allowMM 
            }
        };
    }

    for( size_t i = 0 ; i < threadNum; ++i )
    {
        if( threads[i]->joinable() )
            threads[i]->join();
        delete threads[i];
    }

    std::cerr << ssMsg.rdbuf();

    return 0;
}

std::map< std::string, std::string > get_options( int& argc, char** argv )
{
    std::map< std::string, std::string > args;
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );
    options.add_options()( "help,h" , "Print help messages" )
        ( "fastq,f"    , boost::program_options::value< std::string >()->required()            , "Set input fastq file" )
        ( "index,i"    , boost::program_options::value< std::string >()->required()            , "Set input index file" )
        ( "genome,g"   , boost::program_options::value< std::string >()->default_value( "" )   , "Set input genome fasta file" )
        ( "thread,p"   , boost::program_options::value< std::string >()->default_value( "1" )  , "Set number of threads" )
        ( "minLen,l"   , boost::program_options::value< std::string >()->default_value( "18")  , "Set min length for alignment seed" )
        ( "mismatch,m" , boost::program_options::value< bool >()       ->default_value( false ) , "Set min length for alignment seed" )
        ;
    try
    {
        boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), op );
        if( op.count( "help" ))
        {
            std::cout << "\n" << options << "\n";
            exit(0);
        }
        boost::program_options::notify( op );
    }
    catch( boost::program_options::error& error )
    {
        std::cerr << "\nERROR: " << error.what() << "\n";
        std::cerr << options << "\n";
        exit(1);
    }
    for( auto& arg : op )
        args.emplace( arg.first, arg.first != "mismatch"
                ? arg.second.as< std::string >()
                : arg.second.as< bool >() ? "1" : "0"
        );
    return args;
}
