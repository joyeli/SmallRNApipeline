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

using ArgsType = std::tuple< std::vector< std::string >, std::string, std::string, std::size_t, double, double, double, bool >;
using FileTableType = std::vector< std::vector< std::vector< std::string >>>; // dataset : file : tab
using FileType = std::vector< std::vector< std::string >>; // file : tab

std::vector< std::string > explode( const std::string& in, const std::string& spliter )
{
    std::vector< std::string > content;
    boost::iter_split( content, in, boost::algorithm::first_finder( spliter ));
    return content;
}

std::vector< std::string > get_names( auto& samples )
{
    std::vector< std::string > res, splits;
    for( auto& smp : samples )
    {
        splits = explode( smp, "/" );
        res.emplace_back( splits[0] );
    }
    return res;
}

FileTableType read_files( auto& samples, auto& difmode, auto& diftype )
{
    std::string line;
    FileTableType res;
    std::vector< std::string > temp;

    bool header = true;
    for( auto& smp : samples )
    {
        header = true;
        FileType file;

        std::fstream filein( smp + "/Difference/" + difmode + "Difference_" + diftype + ".tsv", std::ios::in );
        while( std::getline( filein, line ))
        {
            if( header )
            {
                header = false;
                continue;
            }

            temp = explode( line, "\t" );
            file.emplace_back( temp );
        }

        filein.close();
        res.emplace_back( file );
    }

    return res;
}

void filter_via_arm( auto& sample_tables, auto& mexprge, auto& nexprge )
{
    FileTableType res;
    std::vector< std::string > smp_splits;
    std::vector< std::string > arm_splits;
    std::map< std::string, std::size_t > arm_count;

    for( auto& smp : sample_tables )
    {
        FileType file;
        for( auto& line : smp ) 
        {
            if( line[2] == "Y" &&
              ( nexprge != 0 ? nexprge <= std::stod( line[1] ) : true ) &&
              ( mexprge != 0 ? mexprge >= std::stod( line[1] ) : true ) )
            {
                arm_count.clear();
                smp_splits = explode( line[3], ":" );
                arm_splits = explode( line[4], ":" );

                for( auto& p : arm_splits )
                {
                    if( arm_count.find( p ) == arm_count.end() ) arm_count[ p ] = 0;
                    arm_count[ p ]++;
                }

                for( std::size_t i = 0; i < arm_splits.size(); ++i )
                {
                    if( arm_count[ arm_splits[i] ] == 1 )
                        file.emplace_back( std::vector< std::string >{ line[0], smp_splits[i], arm_splits[i] });
                        //                                              anno      sample         arm
                }
            }
        }
        res.emplace_back( file );
    }
    sample_tables = res;
}

void filter_via_len( auto& sample_tables, auto& mexprge, auto& nexprge, auto& filrlen )
{
    FileTableType res;
    for( auto& smp : sample_tables )
    {
        FileType file;
        for( auto& line : smp ) 
        {
            if(( filrlen != 0 ? filrlen <= std::stoi( line[2] ) : true ) &&
               ( nexprge != 0 ? nexprge <= std::stod( line[1] ) : true ) &&
               ( mexprge != 0 ? mexprge >= std::stod( line[1] ) : true ) )
            {
                file.emplace_back( std::vector< std::string >{ line[0], line[2], line[3], line[4] });
                //                                              anno    lendiff  samp1:2  len1:2
            }
        }
        res.emplace_back( file );
    }
    sample_tables = res;
}

void filter_via_lod( auto& sample_tables, auto& mexprge, auto& nexprge, auto& filrchg )
{
    FileTableType res;
    std::vector< std::string > anno_splits;
    std::vector< std::string > samp_splits;

    for( auto& smp : sample_tables )
    {
        FileType file;
        for( auto& line : smp ) 
        {
            if(( filrchg != 0 ? ( std::stod( line[2] ) > 0
                              ? filrchg <= std::stod( line[2] )
                              : filrchg >= std::stod( line[2] )) : true ) &&
               ( nexprge != 0 ? nexprge <= std::stod( line[1] )  : true ) &&
               ( mexprge != 0 ? mexprge >= std::stod( line[1] )  : true ) )
            {
                anno_splits = explode( line[0], "_" );
                samp_splits = explode( line[3], ":" );
                file.emplace_back( std::vector< std::string >{ anno_splits[0], anno_splits[1], line[2], samp_splits[0], samp_splits[1] });
                //                                              anno            seed            change      sample1         sample2
            }
        }
        res.emplace_back( file );
    }
    sample_tables = res;
}

std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting_arm( auto& sample_tables )
{
    std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting;
    for( std::size_t i = 0; i < sample_tables.size(); ++i )
        for( auto& anno : sample_tables[i] ) counting[ anno[1] ][ anno[0] ][ anno[2] ].emplace( i );
    return counting;
}

std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting_len( auto& sample_tables )
{
    std::vector< std::string > smp_splits, len_splits;
    std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting;

    for( std::size_t i = 0; i < sample_tables.size(); ++i )
        for( auto& anno : sample_tables[i] )
        {
            smp_splits = explode( anno[2], ":" );
            len_splits = explode( anno[3], ":" );
            counting[ smp_splits[0] ][ anno[0] ][ len_splits[0] ].emplace( i );
        }
    return counting;
}

std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting_lod( auto& sample_tables )
{
    std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counting;
    for( std::size_t i = 0; i < sample_tables.size(); ++i )
        for( auto& anno : sample_tables[i] ) counting[ anno[3] ][ anno[0] ][ anno[1] ].emplace( i );
    return counting;
}

std::vector< std::pair< std::size_t, std::string >> counting_detail( auto& dataset_name, auto& counts )
{
    std::string string;
    std::set< std::size_t > sets;
    std::vector< std::pair< std::size_t, std::string >> detail_table;

    for( auto& samp : counts ) for( auto& anno : samp.second )
    {
        sets.clear();
        string = samp.first + "\t" + anno.first;

        for( auto& info : anno.second ) for( auto& set : info.second )
        {
            sets.emplace( set );
            string += "\t" + dataset_name[ set ] + ":" + info.first;
        }

        detail_table.emplace_back( std::make_pair( sets.size(), string ));
    }

    return detail_table;
}

std::vector< std::pair< std::size_t, std::string >> counting_detail_limited( auto& dataset_name, auto& counts )
{
    std::string string;
    std::set< std::size_t > sets;
    std::vector< std::pair< std::size_t, std::string >> detail_table;

    for( auto& samp : counts ) for( auto& anno : samp.second ) for( auto& info : anno.second )
    {
        sets.clear();
        string = samp.first + "\t" + anno.first + "_" + info.first;

        for( auto& set : info.second )
        {
            sets.emplace( set );
            string += "\t" + dataset_name[ set ];
        }

        detail_table.emplace_back( std::make_pair( sets.size(), string ));
    }

    return detail_table;
}

void outputing_detail( auto& output, auto& detail_table )
{
    std::sort( detail_table.begin(), detail_table.end(), [] ( const std::pair< std::size_t, std::string >& a, const std::pair< std::size_t, std::string >& b )
    { 
        if( a.first == b.first )
            return a.second > b.second;
        else
            return a.first > b.first;
    });

    for( auto& out : detail_table ) if( out.first > 1 ) output << out.first << "\t" << out.second << "\n";
}

std::map< std::set< std::string >, std::size_t > counting_report( auto& dataset_name, auto& counts )
{
    std::set< std::string > sets;
    std::map< std::set< std::string >, std::size_t > report_table;

    for( auto& samp : counts ) for( auto& anno : samp.second )
    {
        sets.clear();

        for( auto& info : anno.second ) for( auto& set : info.second )
        {
            sets.emplace( dataset_name[ set ]);
        }

        if( report_table.find( sets ) == report_table.end() )
            report_table[ sets ] = 0;

        report_table[ sets ]++;
    }

    return report_table;
}

std::map< std::set< std::string >, std::size_t > counting_report_limited( auto& dataset_name, auto& counts )
{
    std::set< std::string > sets;
    std::map< std::set< std::string >, std::size_t > report_table;

    for( auto& samp : counts ) for( auto& anno : samp.second ) for( auto& info : anno.second )
    {
        sets.clear();

        for( auto& set : info.second )
        {
            sets.emplace( dataset_name[ set ]);
        }

        if( report_table.find( sets ) == report_table.end() )
            report_table[ sets ] = 0;

        report_table[ sets ]++;
    }

    return report_table;
}

void outputing_report( auto& output, auto& report_table )
{
    std::string samples;
    std::vector< std::pair< std::string, std::size_t >> out_vec;

    for( auto& report : report_table )
    {
        samples = "";
        for( auto& sample : report.first ) samples += "\t" + sample;
        out_vec.emplace_back( std::make_pair( samples, report.second ));
    }

    std::sort( out_vec.begin(), out_vec.end(), [] ( const std::pair< std::string, std::size_t >& a, const std::pair< std::string, std::size_t >& b )
    { 
        if( a.first.length() == b.first.length() )
            return a.second > b.second;
        else
            return a.first.length() > b.first.length();
    });

    for( auto& out : out_vec ) output << out.second << "\t" << out.first << "\n";
}

ArgsType get_options( int& argc, char** argv )
{
    ArgsType args;
    bool limited = false;
    boost::program_options::variables_map op;
    boost::program_options::options_description options( "Options" );
    options.add_options()( "help,h" , "Print help messages" )
        ( "samples,s", boost::program_options::value< std::vector< std::string >>()->multitoken(), "Set paths of dataset with biotype; like: GSE45506/miRNA" )
        ( "difmode,d", boost::program_options::value< std::string >()->required(), "Set mode of difference comparision; like: Arm, Length, Loading" )
        ( "diftype,t", boost::program_options::value< std::string >()->required(), "Set type of difference comparision; like: GMPM, GM, PM, Tailing" )
        ( "filrlen,l", boost::program_options::value< std::size_t >()->default_value( 0 ), "Set length of difference for length compare" )
        ( "filrchg,c", boost::program_options::value< double >()->default_value( 0.0 ), "Set %change of difference for loading compare" )
        ( "mexprge,m", boost::program_options::value< double >()->default_value( 0.0 ), "Set max expression for range filter" )
        ( "nexprge,n", boost::program_options::value< double >()->default_value( 0.0 ), "Set max expression for range filter" )
        ( "limited,f", boost::program_options::bool_switch( &limited ), "Set flag of summary limitation of annotation and information, Off: Mir-1 / On: Mir-1_AAAAAAA" )
        ;
    try {
        boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), op );
        if( op.count( "help" )){ std::cout << "\n" << options << "\n"; exit(0); }
        boost::program_options::notify( op );
    }
    catch( boost::program_options::error& error ){ std::cerr << "\nERROR: " << error.what() << "\n" << options << "\n"; exit(1); }
    for( auto& arg : op ) {
        if( arg.first == "samples" ) std::get<0>( args ) = arg.second.as< std::vector< std::string >>();
        if( arg.first == "difmode" ) std::get<1>( args ) = arg.second.as< std::string >();
        if( arg.first == "diftype" ) std::get<2>( args ) = arg.second.as< std::string >();
        if( arg.first == "filrlen" ) std::get<3>( args ) = arg.second.as< std::size_t >();
        if( arg.first == "filrchg" ) std::get<4>( args ) = arg.second.as< double >();
        if( arg.first == "mexprge" ) std::get<5>( args ) = arg.second.as< double >();
        if( arg.first == "nexprge" ) std::get<6>( args ) = arg.second.as< double >();
    }
    std::get<7>( args ) = limited;
    return args;
}

int main( int argc, char** argv )
{
    ArgsType args( get_options( argc, argv ));

    auto& samples = std::get<0>( args );
    auto& difmode = std::get<1>( args );
    auto& diftype = std::get<2>( args );
    auto& filrlen = std::get<3>( args );
    auto& filrchg = std::get<4>( args );
    auto& mexprge = std::get<5>( args );
    auto& nexprge = std::get<6>( args );
    auto& limited = std::get<7>( args );

    std::vector< std::string > dataset_name = get_names( samples );
    FileTableType sample_tables = read_files( samples, difmode, diftype );

    if( difmode == "Arm"     ) filter_via_arm( sample_tables, mexprge, nexprge );
    if( difmode == "Length"  ) filter_via_len( sample_tables, mexprge, nexprge, filrlen );
    if( difmode == "Loading" ) filter_via_lod( sample_tables, mexprge, nexprge, filrchg );

    //          sample                  anno                    info                dataset
    std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>> counts =
    (
        difmode == "Arm"     ? counting_arm( sample_tables ) :
        difmode == "Length"  ? counting_len( sample_tables ) :
        difmode == "Loading" ? counting_lod( sample_tables ) :
        std::map< std::string, std::map< std::string, std::map< std::string, std::set< std::size_t >>>>()
    );

    std::ofstream output;
    std::string output_name = "_" + difmode + "_" + diftype
             + ( filrlen == 0 ? "" : "_len-" + std::to_string( filrlen ))
             + ( filrchg == 0 ? "" : "_chg-" + std::to_string( filrchg ))
             + ( mexprge == 0 && nexprge == 0 ? "" :
                 mexprge == 0 && nexprge != 0 ? "_-" + std::to_string( std::size_t( nexprge )) :
                 mexprge != 0 && nexprge == 0 ? "_"  + std::to_string( std::size_t( mexprge )) + "-" :
                 "_"  + std::to_string( std::size_t( mexprge )) + "-" + std::to_string( std::size_t( nexprge )))
             + ( limited ? "_limited" : "" )
             + ".txt";

    std::vector< std::pair< std::size_t, std::string >> detail_table =
    (
        limited
            ? counting_detail_limited( dataset_name, counts )
            : counting_detail( dataset_name, counts )
    );

    output.open( "CountsDetail" + output_name );
    outputing_detail( output, detail_table );
    output.close();

    std::map< std::set< std::string >, std::size_t > report_table =
    (
        limited
            ? counting_report_limited( dataset_name, counts )
            : counting_report( dataset_name, counts )
    );

    output.open( "CountsReport" + output_name );
    outputing_report( output, report_table );
    output.close();

    return 0;
}
