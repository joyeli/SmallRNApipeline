#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <set>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

int main( int argc, char** argv )
{
    bool noarm;
    std::size_t sampe;
    std::string input, line;
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );
    options.add_options()( "help,h" , "Print help messages" )
        ( "input,i", boost::program_options::value< std::string >( &input )->required(), "Set input tsv file under \"ValPlot\" folder" )
        ( "sampe,s", boost::program_options::value< std::size_t >( &sampe )->required(), "Set sample column#, 0 as all samples but not in the rank" )
        ( "noarm,a", boost::program_options::bool_switch( &noarm )->default_value( false ), "Set arm merge in the annotation" )
        ;
    try{ boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), op );
        if( op.count( "help" )) { std::cout << "\n" << options << "\n"; exit(0); } boost::program_options::notify( op ); }
    catch( boost::program_options::error& error ){ std::cerr << "\nERROR: " << error.what() << "\n" << options << "\n"; exit(1); }

    std::vector< std::string > annons;
    std::vector< std::string > splits;

    std::map< std::string, std::vector< double >> rankmap;
    std::vector< std::pair< double, std::string >> rankvec;

    std::fstream infile( input, std::ios::in );

    while( std::getline( infile, line ))
    {
        splits = explode( line, "\t" );
        if( splits[0] == "Annotation" ) continue;
        annons = explode( splits[0], "_" );

        if( noarm ) annons[0]
            =( annons[0].substr( annons[0].length() -2, 2 ) == "5p"
            || annons[0].substr( annons[0].length() -2, 2 ) == "3p"
            ?  annons[0].substr( 0, annons[0].length() -3 ) : annons[0] );

        if( rankmap.find( annons[0] ) == rankmap.end() ) rankmap[ annons[0] ] = std::vector< double >( splits.size() -1, 0.0 );
        for( std::size_t i = 1; i < splits.size(); ++i ) rankmap[ annons[0] ][ i-1 ] += std::stod( splits[i] );
    }

    infile.close();

    if( sampe == 0 )
    {
        for( auto& anno : rankmap )
        {
            std::cout << anno.first;
            for( auto& smp : anno.second ) std::cout << "\t" << smp;
            std::cout << "\n";
        }
    }
    else
    {
        for( auto& anno : rankmap ) rankvec.emplace_back( std::make_pair( anno.second[ sampe -1 ], anno.first ));

        std::sort( rankvec.begin(), rankvec.end(), []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
        {
            if( a.first == b.first ) return a.second > b.second;
            else return a.first > b.first;
        });

        for( auto& anno : rankvec ) std::cout << anno.second << "\t" << anno.first << "\n";
    }

    return 0;
}
