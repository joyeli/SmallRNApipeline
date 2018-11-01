#pragma once

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerBoxPlotTest
{

  public:

    GeneTypeAnalyzerBoxPlotTest()
    {}

    static void make_file_link( const std::string& p )
    {
        std::string filename;
        std::vector< std::string > list;
        boost::filesystem::path path( p + "../TailDot/" );

        for( auto& file : boost::filesystem::directory_iterator( path ))
        {
            filename = file.path().filename().string();
            if( filename.substr( filename.length() -3  ) == "tsv" &&
              ( filename.length() < 11 ||(
              ( filename.length() > 11 && filename.substr( 0, 14 ) != "Heterorgeneity" ) && 
              ( filename.length() > 11 && filename.substr( filename.length() -11 ) != "isomiRs.tsv" )) ))
                list.emplace_back( filename );
        }

        for( auto& file : list )
            if( !boost::filesystem::exists( p + file ))
                 boost::filesystem::create_symlink(( "../TailDot/" + file ).c_str(), ( p + file ).c_str() );
    }


    static void output_boxplot_visualization( const std::string& output_name )
    {
        make_file_link( output_name );

        std::ofstream output( output_name + "index.php" );
        output.close();
    }

};

} // end of namespace algorithm
} // end of namespace ago
