#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <boost/program_options.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <pokemon/format/annotation_raw_bed.hpp>

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

std::map< std::string, std::string > get_options( int& argc, char** argv )
{
    std::map< std::string, std::string > args;
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );
    options.add_options()( "help,h" , "Print help messages" )
        ( "gnome,g", boost::program_options::value< std::string >()->required(),"Set human genome fasta file" )
        ( "arche,a", boost::program_options::value< std::string >()->required(),"Set archive input file" )
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
    
    auto& gnome = args[ "gnome" ];
    auto& arche = args[ "arche" ];

    std::string line, temp;
    std::map< std::string, std::string > genome_map;
    std::vector< AnnotationRawBed<> > annotation_rawbeds;

    std::ofstream output;
    std::ifstream filein;
    filein.open( gnome );

    while( std::getline( filein, line ))
    {
        if( line.at(0) == '>' )
        {
            temp = line.substr(1);
            genome_map[ temp ] = "";
            continue;
        }

        genome_map[ temp ] += line;
    }

    filein.close();
    filein.open( arche );

    boost::archive::binary_iarchive archein( filein );
    archein & annotation_rawbeds;

    filein.close();
    output.open( arche.substr( 0, arche.size() -4 ) + "_annobed.tsv" );

    output << "Chr\tStart\tEnd\tStrand\tAlignCounts\tReadCounts\tLength\tTailLen\tSeq\tTail\tType\tAnno\tSeed\n";

    for( auto& anno : annotation_rawbeds )
    {
        for( auto& info : anno.annotation_info_ )
        {
            for( int i = 0; i < info.size(); i+=2 )
            {
                output
                    << anno.chromosome_ << "\t"
                    << anno.start_ << "\t"
                    << anno.end_ << "\t"
                    << anno.strand_ << "\t"
                    << anno.multiple_alignment_site_count_ << "\t"
                    << anno.reads_count_ << "\t"
                    << (int)anno.length_ - (int)anno.tail_length_ << "\t"
                    << (int)anno.tail_length_ << "\t"
                    << anno.getReadSeq( genome_map ) << "\t"
                    << anno.getTail() << "\t"
                    << info[i] << "\t"
                    << info[ i+1 ] << "\t"
                    << anno.getReadSeq( genome_map ).substr(1,7) << "\n"
                    ;
            }
        }
    }

    output.close();
    return 0;
}
