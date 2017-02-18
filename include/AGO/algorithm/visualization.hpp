#pragma once
#include <boost/filesystem.hpp>

namespace ago {
namespace algorithm {

class Visualization
{
  public:

    void make_biotype( 
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
	{
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string biotype_path( output_path + "/Visualization/BioType" );
        boost::filesystem::create_directory( biotype_path );

        std::ofstream output;

        for( auto& result_type : all_results )
        {
            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".LenDist" )
            {
                if( gmpm[2] == "read" )
                {
                    int flag = 0;
                    output.open( biotype_path + "/All_Sample_" + gmpm[1] + "_readcount.tsv" );

                    for( std::map< std::string, std::vector< double >>::iterator biotype = result_type.second.begin()->second.begin();
                            biotype != result_type.second.begin()->second.end(); ++biotype )
                    {
                        if( biotype->first == "SUM_ANNO" )
                            continue;

                        if( flag == 0 )
                            output << "BioType";
                        else
                            output << biotype->first;

                        for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                                sample != result_type.second.end(); ++sample )
                        {
                            std::map< std::string, std::vector< double >>::iterator value = sample->second.find( biotype->first );

                            if( flag == 0 )
                                output << "\t" << sample->first;
                            else
                            {
                                output.precision(0);

                                if( value != sample->second.end() && !value->second.empty() )
                                    output << std::fixed << "\t" << value->second[ value->second.size() -1 ];
                                else
                                    output << "\t0";
                            }
                        }
                        flag++;
                        output << "\n";
                    }
                    output.close();
                }
                else
                {
                    int flag = 0;
                    output.open( biotype_path + "/All_Sample_" + gmpm[1] + "_ppm.tsv" );

                    for( std::map< std::string, std::vector< double >>::iterator biotype = result_type.second.begin()->second.begin();
                            biotype != result_type.second.begin()->second.end(); ++biotype )
                    {
                        if( biotype->first == "SUM_ANNO" )
                            continue;

                        if( flag == 0 )
                            output << "BioType";
                        else
                            output << biotype->first;

                        for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                                sample != result_type.second.end(); ++sample )
                        {
                            std::map< std::string, std::vector< double >>::iterator value = sample->second.find( biotype->first );

                            if( flag == 0 )
                                output << "\t" << sample->first;
                            else
                            {
                                output.precision(0);

                                if( value != sample->second.end() && !value->second.empty() )
                                    output << std::fixed << "\t" << value->second[ value->second.size() -1 ];
                                else
                                    output << "\t0";
                            }
                        }

                        flag++;
                        output << "\n";
                    }
                    output.close();
                }
            }
        }

        for( auto& result_type : all_results )
        {
            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".LenDist" )
            {
                if( gmpm[2] == "read" )
                {
                    for( auto& result_sample : result_type.second)
                    {
                        int flag = 0;
                        output.open( biotype_path + "/" + result_sample.first + "_" + gmpm[1] + "_readcount.tsv" );

                        for( auto& result_table : result_sample.second )
                        {
                            if( result_table.first == "SUM_ANNO" )
                                continue;

                            if( flag == 0 )
                                output << "BioType";
                            else
                                output << result_table.first;

                            for( int length = 0; length < result_table.second.size(); ++length )
                            {
                                if( length == result_table.second.size()-1 )
                                {}
                                else
                                {
                                    output.precision(0);
                                    output << std::fixed << "\t" << result_table.second[ length ];
                                }
                            }

                            flag++;
                            output << "\n";
                        }

                        output.close();
                    }
                }
                else
                {
                    for( auto& result_sample : result_type.second)
                    {
                        int flag = 0;
                        output.open( biotype_path + "/" + result_sample.first + "_" + gmpm[1] + "_ppm.tsv" );

                        for( auto& result_table : result_sample.second )
                        {
                            if( result_table.first == "SUM_ANNO" )
                                continue;

                            if( flag == 0 )
                                output << "BioType";
                            else
                                output << result_table.first;

                            for( int length = 0; length < result_table.second.size(); ++length )
                            {
                                if( length == result_table.second.size()-1 )
                                {}
                                else
                                {
                                    output.precision(0);
                                    output << std::fixed << "\t" << result_table.second[ length ];
                                }
                            }

                            flag++;
                            output << "\n";
                        }

                        output.close();
                    }
                }
            }
        }
    }

    void make_lendist(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string lendist_path( output_path + "/Visualization/LenDist" );
        boost::filesystem::create_directory( lendist_path );
        std::ofstream output;

        for( auto& result_type : all_results )
        {
            std::map< int, int > length_index;

            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".LenDist" )
            {
                int flag = 0;
                for( auto& biotype : result_type.second.begin()->second )
                {
                    for( auto& result_sample : result_type.second)
                    {
                        if( flag == 0 )
                        {
                            std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                            for( auto& length : anno->second )
                                length_index.emplace( length, 0 );
                        }
                        flag++;
                    }
                }

                if( gmpm[2] == "read" )
                {
                    output.open( lendist_path + "/All_Sample_" + gmpm[1] + "_readcount.tsv" );
                    output << "miRNA";

                    for( auto& length : length_index )
                    {
                        if( length.first != 0 )
                            output << "\t" << length.first; 
                    }

                    output << "\n";

                    std::vector< std::vector< int >> length_index_check;

                    int it = 0;
                    for( auto& biotype : result_type.second.begin()->second )
                    {
                        if( biotype.first == "SUM_ANNO" )
                            continue;

                        for( auto& result_sample : result_type.second)
                        {
                            std::vector< std::string > anno_check;
                            boost::iter_split( anno_check, biotype.first, boost::algorithm::first_finder( "_" ));

                            if( anno_check[0] == ".LenDist" )
                            {
                                std::vector< int > length_index_check_new;
                                std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                                for( auto& length : anno->second )
                                {
                                    length_index_check_new.push_back( length );
                                }

                                length_index_check.push_back( length_index_check_new );
                                continue;
                            }

                            if( biotype.first != "miRNA" )
                                continue;

                            output << result_sample.first;
                            std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                                if( gmpm[1] == "GMPM" )
                                {
                                    int length = 0;

                                    for( auto& index : length_index )
                                    {
                                        length++;
                                    }


                                    length = 0;

                                    for( auto& index : length_index )
                                    {
                                        length++;
                                    }
                                }

                            int length = 0;
                            for( auto& index : length_index )
                            {
                                if( index.first == 0 )
                                    continue;

                                if( index.first != length_index_check[it][length] )
                                {
                                    output << "\t" << 0;
                                    continue;
                                }

                                output.precision(0);
                                output << std::fixed << "\t" << anno->second[ length ];
                                length++;
                            }

                            output << "\n";
                            it++;
                        }
                    }
                    output.close();
                }
                else
                {
                    output.open( lendist_path + "/All_Sample_" + gmpm[1] + "_ppm.tsv" );
                    output << "miRNA";

                    for( auto& length : length_index )
                    {
                        if( length.first != 0 )
                            output << "\t" << length.first; 
                    }

                    output << "\n";

                    std::vector< std::vector< int >> length_index_check;

                    int it = 0;
                    for( auto& biotype : result_type.second.begin()->second )
                    {
                        if( biotype.first == "SUM_ANNO" )
                            continue;

                        for( auto& result_sample : result_type.second)
                        {
                            std::vector< std::string > anno_check;
                            boost::iter_split( anno_check, biotype.first, boost::algorithm::first_finder( "_" ));

                            if( anno_check[0] == ".LenDist" )
                            {
                                std::vector< int > length_index_check_new;
                                std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                                for( auto& length : anno->second )
                                    length_index_check_new.push_back( length );

                                length_index_check.push_back( length_index_check_new );
                                continue;
                            }

                            if( biotype.first != "miRNA" )
                                continue;

                            output << result_sample.first;
                            std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                            int length = 0;
                            for( auto& index : length_index )
                            {
                                if( index.first == 0 )
                                    continue;

                                if( index.first != length_index_check[it][length] )
                                {
                                    output << "\t" << 0;
                                    continue;
                                }

                                output.precision(0);
                                output << std::fixed << "\t" << anno->second[ length ];
                                length++;
                            }

                            output << "\n";
                            it++;
                        }
                    }
                    output.close();
                }
            }
        }
    }

    void make_mirdist(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string mirdist_path( output_path + "/Visualization/MirDist" );
        boost::filesystem::create_directory( mirdist_path );
        std::ofstream output;

        for( auto& result_type : all_results )
        {
            std::map< int, int > length_index;
            std::map< std::string, double > map_index;
            std::vector< std::pair< std::string, double >> mir_index;

            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".MirDist" )
            {
                if( gmpm[2] == "ppm" )
                {
                    for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                            sample != result_type.second.end(); ++sample )
                    {
                        for( std::map< std::string, std::vector< double >>::iterator biotype = sample->second.begin();
                                biotype != sample->second.end(); ++biotype )
                        {
                            map_index.emplace( biotype->first, 0 );
                        }
                    }

                    for( auto& anno : map_index )
                    {
                        if( anno.first == "SUM_ANNO_SUM_ANNO" || anno.first == ".MirDist_"+gmpm[1]+"_ppm" )
                            continue;

                        double mir_sum(0);

                        for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                                sample != result_type.second.end(); ++sample )
                        {
                            std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );

                            if( value != sample->second.end() && !value->second.empty() )
                            {
                                mir_sum = value->second[ value->second.size() -1 ] + mir_sum;
                            }
                        }

                        mir_index.push_back( std::make_pair( anno.first, mir_sum ));
                    }

                    std::sort( mir_index.begin(), mir_index.end(), 
                            []( const std::pair< std::string , double >& a, const std::pair< std::string, double >& b )
                            { return a.second > b.second; });

                    int flag = 0;
                    for( auto& biotype : result_type.second.begin()->second )
                    {
                        for( auto& result_sample : result_type.second)
                        {
                            if( flag == 0 )
                            {
                                std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                                for( auto& length : anno->second )
                                {
                                    if( length == 0 )
                                        continue;

                                    length_index.emplace( length, 0 );
                                }
                            }
                            flag++;
                        }
                    }

                    for( auto& length : length_index )
                    {
                        output.open( mirdist_path + "/miRNA_" + gmpm[1] + "_length_" + std::to_string( length.first ) + ".tsv" );

                        int flag = 0;

                        for( auto& anno : mir_index )
                        {
                            if( flag == 0 )
                                output << "miRNA";
                            else
                                output << anno.first;

                            for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                                    sample != result_type.second.end(); ++sample )
                            {
                                std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );

                                int length_it = 0;
                                for( auto& len : sample->second.begin()->second )
                                {
                                    if( len == length.first )
                                        break;

                                    length_it++;

                                    if( length_it == sample->second.begin()->second.size() )
                                        length_it = -1;
                                }

                                if( flag == 0 )
                                    output << "\t" << sample->first;
                                else
                                {
                                    if( value != sample->second.end() && !value->second.empty() )
                                    {
                                        if( length_it != -1 )
                                        {
                                            output.precision(0);
                                            output << std::fixed << "\t" << value->second[ length_it ];
                                        }
                                        else
                                        {
                                            output << "\t0";
                                        }
                                    }
                                    else
                                    {
                                        output << "\t0";
                                    }
                                }
                            }

                            if( flag == 0 )
                            {
                                output << "\n" << anno.first;

                                for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                                        sample != result_type.second.end(); ++sample )
                                {
                                    std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );

                                    int length_it = 0;
                                    for( auto& len : sample->second.begin()->second )
                                    {
                                        if( len == length.first )
                                            break;

                                        length_it++;

                                        if( length_it == sample->second.begin()->second.size() )
                                            length_it = -1;
                                    }

                                    if( value != sample->second.end() && !value->second.empty() )
                                    {
                                        if( length_it != -1 )
                                        {
                                            output.precision(0);
                                            output << std::fixed << "\t" << value->second[ length_it ];
                                        }
                                        else
                                        {
                                            output << "\t0";
                                        }
                                    }
                                    else
                                    {
                                        output << "\t0";
                                    }
                                }
                            }

                            flag++;
                            output << "\n";
                        }

                        output.close();
                    }

                }
            }
        }
    }

    void make_dotplot(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string dotplot_path( output_path + "/Visualization/DotPlot" );
        boost::filesystem::create_directory( dotplot_path );
        std::ofstream output;

        std::map< int, int > length_index;
        std::vector< std::string > gmpm_index;
        std::map< std::string, int > sample_index;
        std::map< std::string, double > mir_index;

        for( auto& result_type : all_results )
        {
            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            int check = 0;
            for( auto& gmpm_check : gmpm_index )
                if( gmpm_check == gmpm[1] || gmpm[0] == ".GMPMTR" )
                    check = 1;

            if( check != 0 )
                continue;
            gmpm_index.push_back( gmpm[1] );
        }

        for( auto& result_type : all_results )
        {
            for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                    sample != result_type.second.end(); ++sample )
            {
                sample_index.emplace( sample->first, 0 );
            }

            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".MirDist" && gmpm[1] == "GMPM" && gmpm[2] == "ppm" )
            {
                for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                        sample != result_type.second.end(); ++sample )
                {
                    for( std::map< std::string, std::vector< double >>::iterator biotype = sample->second.begin();
                            biotype != sample->second.end(); ++biotype )
                    {
                        mir_index.emplace( biotype->first, 0 );
                    }
                }

                int flag_check = 0;
                for( auto& biotype : result_type.second.begin()->second )
                {
                    for( auto& result_sample : result_type.second)
                    {
                        if( flag_check == 0 )
                        {
                            std::map< std::string, std::vector< double >>::iterator anno = result_sample.second.find( biotype.first );

                            for( auto& length : anno->second )
                                length_index.emplace( length, 0 );
                        }
                        flag_check++;
                    }
                }
            }
        }

        //	All output but mark abundant miRNA
        for( auto& SAMPLE : sample_index )
        {
            output.open( dotplot_path + "/" + SAMPLE.first + ".tsv" );
            output << SAMPLE.first;

            for( auto& GMPM : gmpm_index )
                output << "\t" << GMPM;

            output << "\n";

            std::map< std::string, std::vector< double >> Check_Abundant;

            for( auto& anno : mir_index )
            {
                if( anno.first == "SUM_ANNO_SUM_ANNO" || anno.first == ".MirDist_GMPM_ppm" )
                    continue;

                std::vector< double > len_vector;

                for( auto& result_type : all_results )
                {
                    std::vector< std::string > gmpm;
                    boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

                    for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                            sample != result_type.second.end(); ++sample )
                    {
                        if( SAMPLE.first == sample->first && gmpm[0] == ".MirDist" && (gmpm[2] == "ppm" || gmpm[2] == "Ratio" ))
                        {
                            std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );
                            int length_it = 0;

                            for( auto& len : sample->second.begin()->second )
                            {
                                if( len == length_index[0] )
                                    break;

                                length_it++;

                                if( length_it == sample->second.begin()->second.size() )
                                    length_it = -1;
                            }

                            if( value != sample->second.end() && !value->second.empty() )
                            {
                                if( length_it != -1 )
                                    len_vector.push_back( value->second[ length_it ] );
                                else
                                    len_vector.push_back( 0 );
                            }
                            else
                                len_vector.push_back( 0 );
                        }

                    }
                }
                Check_Abundant.emplace( anno.first, len_vector );
            }

            int check_flag(0);
            std::string anno_checker;
            std::string anno_seed;
            std::vector< std::pair< double, std::string >> abundant_sorted;

            for( auto& anno_pair : Check_Abundant )
            {
                std::vector< std::string > anno;
                boost::iter_split( anno, anno_pair.first, boost::algorithm::first_finder( "_" ));

                if( check_flag == 0 )
                    anno_checker = anno[0];

                if( anno_checker == anno[0] )
                {
                    abundant_sorted.push_back( std::make_pair( anno_pair.second[0], anno[1] ));
                }
                else
                {
                    std::sort( abundant_sorted.begin(), abundant_sorted.end(), 
                            []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
                            { return a.first > b.first; });

                    for( auto& miRNA : abundant_sorted )
                    {
                        if( miRNA.second != abundant_sorted.begin()->second )
                        {
                            anno_seed = anno_checker + "_" + miRNA.second;
                            output << anno_seed;

                            for( auto& len : Check_Abundant.find( anno_seed )->second )
                            {
                                output.precision(0);
                                output << std::fixed << "\t" << len;
                            }

                            output << "\n";
                        }
                        else
                        {
                            anno_seed = anno_checker + "_" + miRNA.second;
                            output << "*" << anno_seed;

                            for( auto& len : Check_Abundant.find( anno_seed )->second )
                            {
                                output.precision(0);
                                output << std::fixed << "\t" << len;
                            }

                            output << "\n";
                        }
                    }

                    anno_checker = anno[0];
                    abundant_sorted.clear();
                    abundant_sorted.push_back( std::make_pair( anno_pair.second[0], anno[1] ));
                }

                ++check_flag;
            }

            std::sort( abundant_sorted.begin(), abundant_sorted.end(), 
                    []( const std::pair< double, std::string >& a, const std::pair< double, std::string >& b )
                    { return a.first > b.first; });

            for( auto& miRNA : abundant_sorted )
            {
                if( miRNA.second != abundant_sorted.begin()->second )
                {
                    anno_seed = anno_checker + "_" + miRNA.second;
                    output << anno_seed;

                    for( auto& len : Check_Abundant.find( anno_seed )->second )
                    {
                        output.precision(0);
                        output << std::fixed << "\t" << len;
                    }

                    output << "\n";
                }
                else
                {
                    anno_seed = anno_checker + "_" + miRNA.second;
                    output << "*" << anno_seed;

                    for( auto& len : Check_Abundant.find( anno_seed )->second )
                    {
                        output.precision(0);
                        output << std::fixed << "\t" << len;
                    }

                    output << "\n";
                }
            }

            output.close();
        }
    }

    void make_valplot(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string valplot_path( output_path + "/Visualization/ValPlot" );
        boost::filesystem::create_directory( valplot_path );
        std::ofstream output;

        std::map< std::string, int > gmpm_index;
        std::vector< std::string > sample_index;
        std::map< std::string, double > mir_index;

        for( auto& result_type : all_results )
        {
            for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                    sample != result_type.second.end(); ++sample )
            {
                int check = 0;
                for( auto& sample_check : sample_index )
                    if( sample_check == sample->first )
                        check = 1;

                if( check != 0 )
                    continue;
                sample_index.push_back( sample->first );
            }

            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] != ".GMPMTR" && gmpm[0] != ".MirTail" )
                gmpm_index.emplace( gmpm[1], 0 );

            if( gmpm[0] == ".MirDist" && gmpm[1] == "GMPM" && gmpm[2] == "ppm" )
            {
                for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                        sample != result_type.second.end(); ++sample )
                {
                    for( std::map< std::string, std::vector< double >>::iterator biotype = sample->second.begin();
                            biotype != sample->second.end(); ++biotype )
                    {
                        mir_index.emplace( biotype->first, 0 );
                    }
                }
            }
        }

        std::map< std::string, double > all_miRNA;

        for( auto& anno : mir_index )
        {
            if( anno.first == "SUM_ANNO_SUM_ANNO" || anno.first == ".MirDist_GMPM_ppm" )
                continue;

            double sum(0);

            for( auto& result_type : all_results )
            {
                std::vector< std::string > gmpm;
                boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));
                for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                        sample != result_type.second.end(); ++sample )
                {
                    if( gmpm[1] == "GMPM" && gmpm[0] == ".MirDist" && (gmpm[2] == "ppm" || gmpm[2] == "Ratio" ))
                    {
                        std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );

                        if( value != sample->second.end() && !value->second.empty() )
                            sum += value->second[ value->second.size()-1 ];
                    }
                }
            }
            all_miRNA.emplace( anno.first, sum );
        }

        std::string anno_check;
        std::map< std:: string, int > most_abundant;
        std::vector< std::pair< std::string, double >> seed_check;

        for( auto& anno : all_miRNA )
        {
            std::vector< std::string > anno_seed;
            boost::iter_split( anno_seed, anno.first, boost::algorithm::first_finder( "_" ));

            if( anno.first == all_miRNA.begin()->first )
                anno_check = anno_seed[0];

            if( anno_seed[0] != anno_check )
            {
                std::sort( seed_check.begin(), seed_check.end(), 
                        []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
                        { return a.second > b.second; });

                most_abundant.emplace( anno_check + "_" + seed_check[0].first, 0 );

                seed_check.clear();
                anno_check = anno_seed[0];
            }

            seed_check.push_back( std::make_pair( anno_seed[1], anno.second ));
        }

        std::sort( seed_check.begin(), seed_check.end(), 
                []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
                { return a.second > b.second; });

        most_abundant.emplace( anno_check + "_" + seed_check[0].first, 0 );

        for( auto& GMPM : gmpm_index )
        {
            output.open( valplot_path + "/" + GMPM.first + ".tsv" );
            output << "miRNA";

            for( auto& SAMPLE : sample_index )
                output << "\t" << SAMPLE;

            output << "\n";

            for( auto& anno : mir_index )
            {
                if( anno.first == "SUM_ANNO_SUM_ANNO" || anno.first == ".MirDist_GMPM_ppm" )
                    continue;

                if( most_abundant.find( anno.first ) != most_abundant.end() )
                    output << "*" << anno.first;
                else
                    output << anno.first;

                for( auto& result_type : all_results )
                {
                    std::vector< std::string > gmpm;
                    boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

                    for( std::map< std::string, std::map< std::string, std::vector< double >>>::iterator sample = result_type.second.begin();
                            sample != result_type.second.end(); ++sample )
                    {
                        if( GMPM.first == gmpm[1] && gmpm[0] == ".MirDist" && (gmpm[2] == "ppm" || gmpm[2] == "Ratio" ))
                        {
                            std::map< std::string, std::vector< double >>::iterator value = sample->second.find( anno.first );

                            if( value != sample->second.end() && !value->second.empty() )
                            {
                                output.precision(0);
                                output << std::fixed << "\t" << value->second[ value->second.size()-1 ];
                            }
                            else
                                output << "\t0";
                        }
                    }
                }
                output << "\n";
            }
            output.close();
        }
    }

    void make_lenplus(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string lenplus_path( output_path + "/Visualization/LenPlus" );
        boost::filesystem::create_directory( lenplus_path );
        std::ofstream output;

        std::vector< std::string > index({ "A_Tail", "C_Tail", "G_Tail", "T_Tail", "Other_Tail", "GM" });
        std::vector< std::string > header;

        for( auto& result_type : all_results )
        {
            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".GMPMTR" && gmpm[1] == "ppm" )
                for( auto& length : result_type.second.begin()->second )
                    if( length.first != ".GMPMTR_ppm" )
                        header.push_back( length.first );
        }

        for( auto& result_type : all_results )
        {
            std::vector< std::string > gmpm;
            boost::iter_split( gmpm, result_type.first, boost::algorithm::first_finder( "_" ));

            if( gmpm[0] == ".GMPMTR" )
            {
                if( gmpm[1] == "ppm" )
                {
                    for( auto& result_sample : result_type.second)
                    {
                        output.open( lenplus_path + "/" + result_sample.first + "_ppm" + ".tsv" );
                        output << "Length";

                        for( auto& length : header )
                            output << "\t" << length;

                        output << "\n" << index[ index.size()-1 ];

                        for( auto& length : result_sample.second )
                        {
                            if( length.first == ".GMPMTR_ppm" )
                                continue;

                            output.precision(0);
                            output << std::fixed << "\t" << length.second[ index.size()-1 ];
                        }

                        output << "\n";						

                        for( int i = 0; i < index.size()-1; i++ )
                        {
                            output << index[i];

                            int flag(0);
                            for( auto& length : result_sample.second )
                            {
                                if( flag == 0 )
                                {
                                    ++flag;
                                    continue;
                                }

                                ++flag;
                                output.precision(0);
                                output << std::fixed << "\t" << length.second[i];
                            }
                            output << "\n";
                        }
                        output.close();
                    }
                }
                else
                {
                    for( auto& result_sample : result_type.second)
                    {
                        output.open( lenplus_path + "/" + result_sample.first + "_readcount" + ".tsv" );
                        output << "Length";

                        for( auto& length : header )
                            output << "\t" << length;

                        output << "\n" << index[ index.size()-1 ];

                        for( auto& length : result_sample.second )
                        {
                            if( length.first == ".GMPMTR_read_count" )
                                continue;

                            output.precision(0);
                            output << std::fixed << "\t" << length.second[ index.size()-1 ];
                        }

                        output << "\n";						

                        for( int i = 0; i < index.size()-1; i++ )
                        {
                            output << index[i];

                            int flag(0);
                            for( auto& length : result_sample.second )
                            {
                                if( flag == 0 )
                                {
                                    ++flag;
                                    continue;
                                }

                                ++flag;
                                output.precision(0);
                                output << std::fixed << "\t" << length.second[i];
                            }
                            output << "\n";
                        }
                        output.close();
                    }
                }
            }
        }
    }

    void make_mirtail(
          std::vector< std::pair< std::string, ago::engine::DataPool::anno_len_samples >>& all_results
        , const std::string& output_path
    )
    {
        boost::filesystem::create_directory( output_path + "/Visualization" );
        std::string mirtail_path( output_path + "/Visualization/MirTail" );
        boost::filesystem::create_directory( mirtail_path );
        std::ofstream output;

        std::vector< std::string > index({ "A_Tail", "C_Tail", "G_Tail", "T_Tail", "Other_Tail", "GM" });

        for( auto& result_type : all_results )
        {
            if( result_type.first == ".MirTail_ppm" )
            {
                for( auto& result_sample : result_type.second)
                {
                    std::map< std::string, double > all_miRNA;

                    for( auto& anno : result_sample.second )
                    {
                        if( anno.first != ".MirTail_ppm" )
                        {
                            std::vector< std::string > anno_len;
                            boost::iter_split( anno_len, anno.first, boost::algorithm::first_finder( ":" ));

                            double value(0);

                            for( auto& tail : anno.second )
                                value += tail;

                            if( all_miRNA.find( anno_len[0] ) != all_miRNA.end() )
                                all_miRNA.find( anno_len[0] )->second += value;
                            else
                                all_miRNA.emplace( anno_len[0], value );
                        }
                    }

                    std::string anno_check;
                    std::map< std:: string, int > most_abundant;
                    std::vector< std::pair< std::string, double >> seed_check;

                    for( auto& anno : all_miRNA )
                    {
                        std::vector< std::string > anno_seed;
                        boost::iter_split( anno_seed, anno.first, boost::algorithm::first_finder( "_" ));

                        if( anno.first == all_miRNA.begin()->first )
                            anno_check = anno_seed[0];

                        if( anno_seed[0] != anno_check )
                        {
                            std::sort( seed_check.begin(), seed_check.end(), 
                                    []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
                                    { return a.second > b.second; });

                            most_abundant.emplace( anno_check + "_" + seed_check[0].first, 0 );

                            seed_check.clear();
                            anno_check = anno_seed[0];
                        }

                        seed_check.push_back( std::make_pair( anno_seed[1], anno.second ));
                    }

                    std::sort( seed_check.begin(), seed_check.end(), 
                            []( const std::pair< std::string, double >& a, const std::pair< std::string, double >& b )
                            { return a.second > b.second; });

                    most_abundant.emplace( anno_check + "_" + seed_check[0].first, 0 );

                    output.open( mirtail_path + "/" + result_sample.first + "_ppm" + ".tsv" );

                    for( auto& anno : result_sample.second )
                    {
                        if( anno.first == ".MirTail_ppm" )
                        {
                            output << "miRNA_Length";

                            for( auto& tail : index )
                                output << "\t" << tail;

                            output << "\n";
                        }
                        else
                        {
                            std::vector< std::string > anno_len;
                            boost::iter_split( anno_len, anno.first, boost::algorithm::first_finder( ":" ));

                            if( most_abundant.find( anno_len[0] ) != most_abundant.end() )
                                output << anno_len[0] << "*:" << anno_len[1];
                            else
                                output << anno.first;

                            for( auto& value : anno.second )
                                output << "\t" << value;

                            output << "\n";
                        }
                    }

                    output.close();
                }
            }
        }
    }

    void make_html( std::string& output_path )
    {
        std::ofstream output;
        boost::filesystem::create_directory( output_path + "/Visualization" );
        output.open( output_path + "/Visualization/index.php", std::ios::out | std::ios::trunc );

        if( !output.is_open() )
            std::cerr << "\nOutput of index.php is not open\n";
        else
        {
            output << "<? " << "\n" << std::flush;
            output << "	Shell_Exec( 'rm /tmp/*' );" << "\n" << std::flush;
            output << "	$GMPM = $_POST['GMPM'];" << "\n" << std::flush;
            output << "	$FGMPM = $_POST['FGMPM'];" << "\n" << std::flush;
            output << "	$isLog = $_POST['isLog'];" << "\n" << std::flush;
            output << "	$RCPPM = $_POST['RCPPM'];" << "\n" << std::flush;
            output << "	$ForceY = $_POST['ForceY'];" << "\n" << std::flush;
            output << "	$Filter = $_POST['Filter'];" << "\n" << std::flush;
            output << "	$Bar_Pie = $_POST['Bar_Pie'];" << "\n" << std::flush;
            output << "	$MaxHight = $_POST['MaxHight'];" << "\n" << std::flush;
            output << "	$ForceMin = $_POST['ForceMin'];" << "\n" << std::flush;
            output << "	$ForceMax = $_POST['ForceMax'];" << "\n" << std::flush;
            output << "	$TSV_File = $_POST['TSV_File'];" << "\n" << std::flush;
            output << "	$TSV_File2 = $_POST['TSV_File2'];" << "\n" << std::flush;
            output << "	$Top_miRNA = $_POST['Top_miRNA'];" << "\n" << std::flush;
            output << "	$FilterMin = $_POST['FilterMin'];" << "\n" << std::flush;
            output << "	$FilterMax = $_POST['FilterMax'];" << "\n" << std::flush;
            output << "	$Chart_Type = $_POST['Chart_Type'];" << "\n" << std::flush;
            output << "	$isAbundant = $_POST['isAbundant'];" << "\n" << std::flush;
            output << "	$mirDistType = $_POST['mirDistType'];" << "\n" << std::flush;
            output << "	$miRNA_Select = $_POST['miRNA_Select'];" << "\n" << std::flush;
            output << "?>" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "<!DOCTYPE html>" << "\n" << std::flush;
            output << "<html>" << "\n" << std::flush;
            output << "	<meta charset='utf-8'>" << "\n" << std::flush;
            output << "	<body>" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "	<? " << "\n" << std::flush;
            output << "#<!--================== Chart Type ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/d3.min.js></script>';" << "\n" << std::flush;
            output << "		echo '<link href=http://140.113.239.177/joye/ForAgoSorting/lib/svg0331.css rel=stylesheet type=text/css>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		echo '<select name=Chart_Type onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "		echo '<option '; if($Chart_Type=='') echo 'selected'; echo '>Chart Type</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		$Folder = Shell_Exec( 'ls -d */' );" << "\n" << std::flush;
            output << "		$Folder_List = Explode( \"\\n\", $Folder );" << "\n" << std::flush;
            output << "		$Folder_Size = Count( $Folder_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		For( $i = 0; $i < $Folder_Size-1; ++$i )" << "\n" << std::flush;
            output << "		{" << "\n" << std::flush;
            output << "			echo '<option value='.$Folder_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type == $Folder_List[$i] )" << "\n" << std::flush;
            output << "				echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			echo '>'.$Folder_List[$i].'</option>';" << "\n" << std::flush;
            output << "		}" << "\n" << std::flush;
            output << "		echo \"</select>" << "\n" << std::flush;
            output << "			</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== GMPM ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "		if( $Chart_Type != 'ValPlot/' && $Chart_Type != 'LenPlus/' && $Chart_Type != 'MirTail/' )" << "\n" << std::flush;
            output << "		{" << "\n" << std::flush;
            output << "			echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			echo '<select name=GMPM onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "			echo '<option '; if($GMPM=='') echo 'selected'; echo '>GM or PM</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			$GMPM_List = array('GMPM', 'GM', 'PM', 'Tailing');" << "\n" << std::flush;
            output << "			$GMPM_Size = Count( $GMPM_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			For( $i = 0; $i < $GMPM_Size; ++$i )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '<option value='.$GMPM_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $GMPM == $GMPM_List[$i] )" << "\n" << std::flush;
            output << "					echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '>' . $GMPM_List[$i] . '</option>';" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "			echo \"</select>" << "\n" << std::flush;
            output << "				<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='MaxHight' value='$MaxHight' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "				<input type='hidden' name='mirDistType' value='$mirDistType' />" << "\n" << std::flush;
            output << "				</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type != 'MirDist/' && $Chart_Type != 'DotPlot/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.js></script>';" << "\n" << std::flush;
            output << "				echo '<link href=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Single TSV File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep _'.$GMPM.'_ | grep .tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $TSV_File == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $Chart_Type == 'BioType/' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<style>';" << "\n" << std::flush;
            output << "					echo 'text[ style=\"opacity: 0; text-anchor: middle;\" ]{';" << "\n" << std::flush;
            output << "					echo 'opacity: 1 !important';" << "\n" << std::flush;
            output << "					echo '}';" << "\n" << std::flush;
            output << "					echo '</style>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== BioType ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<input type=radio name=Bar_Pie value=Bar onchange=this.form.submit(); ';" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Bar' ) echo 'checked';" << "\n" << std::flush;
            output << "					echo ' >Bar</input>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<input type=radio name=Bar_Pie value=Pie onchange=this.form.submit(); ';" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Pie' ) echo 'checked';" << "\n" << std::flush;
            output << "					echo ' >Pie</input>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \"<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "						<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "						<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "						</form><br/>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$File = $Chart_Type.$TSV_File;" << "\n" << std::flush;
            output << "					$File_Tag = Explode( '_', $TSV_File );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Bar' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						if( $File_Tag[0] == 'All' )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							echo \"<svg id='bar'></svg>" << "\n" << std::flush;
            output << "								<script>" << "\n" << std::flush;
            output << "									d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });" << "\n" << std::flush;
            output << "										var data = new Array();" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										var id = 0;" << "\n" << std::flush;
            output << "										length.forEach( function(d) {" << "\n" << std::flush;
            output << "											var sample = new Object();" << "\n" << std::flush;
            output << "											var value = new Array();" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											sample.key = String( [d] );" << "\n" << std::flush;
            output << "											sample.values = value;" << "\n" << std::flush;
            output << "										" << "\n" << std::flush;
            output << "											data[id] = sample;" << "\n" << std::flush;
            output << "											id++;" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "											id = 0;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "												var value = new Object();" << "\n" << std::flush;
            output << "												value.x = d['BioType'];" << "\n" << std::flush;
            output << "												value.y = Number(d[key]);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "												data[id].values.push( value );" << "\n" << std::flush;
            output << "												id++;" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										nv.addGraph({" << "\n" << std::flush;
            output << "											generate: function() {" << "\n" << std::flush;
            output << "												var width = nv.utils.windowSize().width," << "\n" << std::flush;
            output << "													height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												var chart = nv.models.multiBarChart()" << "\n" << std::flush;
            output << "													.margin({'left':70})" << "\n" << std::flush;
            output << "													.staggerLabels(true)" << "\n" << std::flush;
            output << "													.width(width)" << "\n" << std::flush;
            output << "													.height(height);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												chart.dispatch.on('renderEnd', function(){" << "\n" << std::flush;
            output << "													console.log('Render Complete');" << "\n" << std::flush;
            output << "												});" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "												var svg = d3.select('#bar').datum(data);" << "\n" << std::flush;
            output << "												console.log('calling chart');" << "\n" << std::flush;
            output << "												svg.transition().duration(0).call(chart);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												return chart;" << "\n" << std::flush;
            output << "											}," << "\n" << std::flush;
            output << "											callback: function(graph) {" << "\n" << std::flush;
            output << "												nv.utils.windowResize(function() {" << "\n" << std::flush;
            output << "													var width = nv.utils.windowSize().width;" << "\n" << std::flush;
            output << "													var height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "													graph.width(width).height(height);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "													d3.select('#bar')" << "\n" << std::flush;
            output << "														.attr('width', width)" << "\n" << std::flush;
            output << "														.attr('height', height)" << "\n" << std::flush;
            output << "														.transition().duration(0)" << "\n" << std::flush;
            output << "														.call(graph);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												});" << "\n" << std::flush;
            output << "											}" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "								</script>\";" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "						else" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							echo \"<svg id='bar'></svg>" << "\n" << std::flush;
            output << "								<script>" << "\n" << std::flush;
            output << "									d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });" << "\n" << std::flush;
            output << "										var data = new Array();" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "											d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "												return {" << "\n" << std::flush;
            output << "													x: Number(key)," << "\n" << std::flush;
            output << "													y: +d[key]" << "\n" << std::flush;
            output << "												};" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "											var sample = new Object();" << "\n" << std::flush;
            output << "											sample.key = d['BioType'];" << "\n" << std::flush;
            output << "											sample.values = d.len;" << "\n" << std::flush;
            output << "											data.push( sample );" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										nv.addGraph({" << "\n" << std::flush;
            output << "											generate: function() {" << "\n" << std::flush;
            output << "												var width = nv.utils.windowSize().width," << "\n" << std::flush;
            output << "													height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												var chart = nv.models.multiBarChart()" << "\n" << std::flush;
            output << "													.margin({'left':70})" << "\n" << std::flush;
            output << "													.width(width)" << "\n" << std::flush;
            output << "													.height(height);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												chart.dispatch.on('renderEnd', function(){" << "\n" << std::flush;
            output << "													console.log('Render Complete');" << "\n" << std::flush;
            output << "												});" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "												var svg = d3.select('#bar').datum(data);" << "\n" << std::flush;
            output << "												console.log('calling chart');" << "\n" << std::flush;
            output << "												svg.transition().duration(0).call(chart);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												return chart;" << "\n" << std::flush;
            output << "											}," << "\n" << std::flush;
            output << "											callback: function(graph) {" << "\n" << std::flush;
            output << "												nv.utils.windowResize(function() {" << "\n" << std::flush;
            output << "													var width = nv.utils.windowSize().width;" << "\n" << std::flush;
            output << "													var height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "													graph.width(width).height(height);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "													d3.select('#bar')" << "\n" << std::flush;
            output << "														.attr('width', width)" << "\n" << std::flush;
            output << "														.attr('height', height)" << "\n" << std::flush;
            output << "														.transition().duration(0)" << "\n" << std::flush;
            output << "														.call(graph);" << "\n" << std::flush;
            output << "									" << "\n" << std::flush;
            output << "												});" << "\n" << std::flush;
            output << "											}" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "								</script>\";" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Pie' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						if( $File_Tag[0] == 'All' )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "							$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "							$inFile_Line = Explode( \"\\t\", $inFile_Lines[0] );" << "\n" << std::flush;
            output << "							$idCount = Count( $inFile_Line );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							For( $i = 1; $i < $idCount; ++$i )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								echo '<svg id=pie'.$i.'></svg>';" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							echo \"<script>" << "\n" << std::flush;
            output << "									d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });" << "\n" << std::flush;
            output << "										var data = new Array();" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										var id = 0;" << "\n" << std::flush;
            output << "										length.forEach( function(d) {" << "\n" << std::flush;
            output << "											var sample = new Object();" << "\n" << std::flush;
            output << "											var value = new Array();" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											sample.key = String( [d] );" << "\n" << std::flush;
            output << "											sample.values = value;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											data[id] = sample;" << "\n" << std::flush;
            output << "											id++;" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "											id = 0;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											d.len = length.map( function(c) {" << "\n" << std::flush;
            output << "												var value = new Object();" << "\n" << std::flush;
            output << "												value.key = d['BioType'];" << "\n" << std::flush;
            output << "												value.y = d[c];" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "												data[id].values.push( value );" << "\n" << std::flush;
            output << "												id++;" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										var svg_id = 0;" << "\n" << std::flush;
            output << "										data.forEach( function(d) {" << "\n" << std::flush;
            output << "											console.log( d.key );" << "\n" << std::flush;
            output << "											console.log( d.values );" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											nv.addGraph({" << "\n" << std::flush;
            output << "												generate: function() {" << "\n" << std::flush;
            output << "													svg_id++;" << "\n" << std::flush;
            output << "										" << "\n" << std::flush;
            output << "													var chart = nv.models.pieChart()" << "\n" << std::flush;
            output << "														.x(function(c) { return c.key })" << "\n" << std::flush;
            output << "														.y(function(c) { return c.y })" << "\n" << std::flush;
            output << "														.donut(true)" << "\n" << std::flush;
            output << "														.width('500px')" << "\n" << std::flush;
            output << "														.height('500px')" << "\n" << std::flush;
            output << "														.cornerRadius(10)" << "\n" << std::flush;
            output << "														.showLabels(true);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "													chart.title( d.key );" << "\n" << std::flush;
            output << "													chart.pie.donutLabelsOutside(true).donut(true);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "													var svg = d3.select('#pie' + svg_id).datum( d.values );" << "\n" << std::flush;
            output << "													console.log('calling chart');" << "\n" << std::flush;
            output << "													svg.transition().duration(1200).call(chart);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "													return chart;" << "\n" << std::flush;
            output << "												}" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "								</script>\";" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "						else" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							echo '<br/>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "							$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "							$idCount = Count( $inFile_Lines );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							For( $i = 1; $i < $idCount-1; ++$i )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								echo '<svg id=pie'.$i.'></svg>';" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							echo \"<script>" << "\n" << std::flush;
            output << "									d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "										var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'BioType'; });" << "\n" << std::flush;
            output << "										var svg_id = 0;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "											d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "												return {" << "\n" << std::flush;
            output << "													key: Number(key)," << "\n" << std::flush;
            output << "													y: +d[key]" << "\n" << std::flush;
            output << "												};" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											nv.addGraph({" << "\n" << std::flush;
            output << "												generate: function() {" << "\n" << std::flush;
            output << "													svg_id++;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "													var chart = nv.models.pieChart()" << "\n" << std::flush;
            output << "														.x(function(c) { return c.key })" << "\n" << std::flush;
            output << "														.y(function(c) { return c.y })" << "\n" << std::flush;
            output << "														.donut(true)" << "\n" << std::flush;
            output << "														.width('500px')" << "\n" << std::flush;
            output << "														.height('500px')" << "\n" << std::flush;
            output << "														.cornerRadius(10)" << "\n" << std::flush;
            output << "														.showLabels(true);" << "\n" << std::flush;
            output << "										" << "\n" << std::flush;
            output << "													chart.title( d['BioType'] );" << "\n" << std::flush;
            output << "													chart.pie.donutLabelsOutside(true).donut(true);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "													var svg = d3.select('#pie' + svg_id).datum( d.len );" << "\n" << std::flush;
            output << "													console.log('calling chart');" << "\n" << std::flush;
            output << "													svg.transition().duration(1200).call(chart);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "													return chart;" << "\n" << std::flush;
            output << "												}" << "\n" << std::flush;
            output << "											});" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								</script>\";" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $Chart_Type == 'LenDist/' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== LenDist ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<input type=radio name=Bar_Pie value=Bar onchange=this.form.submit(); ';" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Bar' ) echo 'checked';" << "\n" << std::flush;
            output << "					echo ' >Bar</input>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '<input type=radio name=Bar_Pie value=Pie onchange=this.form.submit(); ';" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Pie' ) echo 'checked';" << "\n" << std::flush;
            output << "					echo ' >Pie</input>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \"" << "\n" << std::flush;
            output << "						<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "						<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "						<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "						</form><br/>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$File = $Chart_Type.$TSV_File;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Bar' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						echo \"<svg id='bar'></svg>" << "\n" << std::flush;
            output << "							<script>" << "\n" << std::flush;
            output << "								d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "								var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA'; });" << "\n" << std::flush;
            output << "								var data = new Array();" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "								tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "									d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "										return {" << "\n" << std::flush;
            output << "											x: Number(key)," << "\n" << std::flush;
            output << "											y: +d[key]" << "\n" << std::flush;
            output << "										};" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "									var sample = new Object();" << "\n" << std::flush;
            output << "									sample.key = d['miRNA'];" << "\n" << std::flush;
            output << "									sample.values = d.len;" << "\n" << std::flush;
            output << "									data.push( sample );" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "								nv.addGraph({" << "\n" << std::flush;
            output << "									generate: function() {" << "\n" << std::flush;
            output << "										var width = nv.utils.windowSize().width," << "\n" << std::flush;
            output << "											height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										var chart = nv.models.multiBarChart()" << "\n" << std::flush;
            output << "											.margin({'left':70})" << "\n" << std::flush;
            output << "											.width(width)" << "\n" << std::flush;
            output << "											.height(height);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										chart.dispatch.on('renderEnd', function(){" << "\n" << std::flush;
            output << "											console.log('Render Complete');" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "										var svg = d3.select('#bar').datum(data);" << "\n" << std::flush;
            output << "										console.log('calling chart');" << "\n" << std::flush;
            output << "										svg.transition().duration(0).call(chart);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										return chart;" << "\n" << std::flush;
            output << "									}," << "\n" << std::flush;
            output << "									callback: function(graph) {" << "\n" << std::flush;
            output << "										nv.utils.windowResize(function() {" << "\n" << std::flush;
            output << "											var width = nv.utils.windowSize().width;" << "\n" << std::flush;
            output << "											var height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "											graph.width(width).height(height);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "											d3.select('#bar')" << "\n" << std::flush;
            output << "												.attr('width', width)" << "\n" << std::flush;
            output << "												.attr('height', height)" << "\n" << std::flush;
            output << "												.transition().duration(0)" << "\n" << std::flush;
            output << "												.call(graph);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "										});" << "\n" << std::flush;
            output << "									}" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "						</script>\";" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $Bar_Pie == 'Pie' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "						$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "						$inFile_Line = Explode( \"\\t\", $inFile_Lines[0] );" << "\n" << std::flush;
            output << "						$idCount = Count( $inFile_Line );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						For( $i = 1; $i < $idCount; ++$i )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							echo '<svg id=pie'.$i.'></svg>';" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						echo \"<script>" << "\n" << std::flush;
            output << "								d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "								var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA'; });" << "\n" << std::flush;
            output << "								var svg_id = 0;" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "								tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "									d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "										return {" << "\n" << std::flush;
            output << "											key: Number(key)," << "\n" << std::flush;
            output << "											y: +d[key]" << "\n" << std::flush;
            output << "										};" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "									nv.addGraph({" << "\n" << std::flush;
            output << "										generate: function() {" << "\n" << std::flush;
            output << "											svg_id++;" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "											var chart = nv.models.pieChart()" << "\n" << std::flush;
            output << "												.x(function(c) { return c.key })" << "\n" << std::flush;
            output << "												.y(function(c) { return c.y })" << "\n" << std::flush;
            output << "												.donut(true)" << "\n" << std::flush;
            output << "												.width('500px')" << "\n" << std::flush;
            output << "												.height('500px')" << "\n" << std::flush;
            output << "												.cornerRadius(10)" << "\n" << std::flush;
            output << "												.showLabels(true);" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "											chart.title( d['miRNA'] );" << "\n" << std::flush;
            output << "											chart.pie.donutLabelsOutside(true).donut(true);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "											var svg = d3.select('#pie' + svg_id).datum( d.len );" << "\n" << std::flush;
            output << "											console.log('calling chart');" << "\n" << std::flush;
            output << "											svg.transition().duration(1200).call(chart);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "											return chart;" << "\n" << std::flush;
            output << "										}" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "						</script>\";" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type == 'MirDist/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '<style>" << "\n" << std::flush;
            output << "						.x.axis path {" << "\n" << std::flush;
            output << "						display: none;" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "					</style>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Top miRNA Select ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=Top_miRNA onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($Top_miRNA=='') echo 'selected'; echo '>Select Top miRNA</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep _'.$GMPM.'_ | grep .tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$Count_File = $Chart_Type . $TSV_List[0];" << "\n" << std::flush;
            output << "				$Count_inFile = File_get_contents( $Count_File );" << "\n" << std::flush;
            output << "				$Count_inFile_Lines = Explode( \"\\n\", $Count_inFile );" << "\n" << std::flush;
            output << "				$Count_miRNA = Count( $Count_inFile_Lines )-2;" << "\n" << std::flush;
            output << "				$Count_miRNA10 = $Count_miRNA % 10;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 10; $i < $Count_miRNA; $i += 10 )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$i.' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $i == $Top_miRNA )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$i.'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<option value='.$Count_miRNA.' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $Count_miRNA == $Top_miRNA )" << "\n" << std::flush;
            output << "					echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '>'.$Count_miRNA.'</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== mirDistType ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=mirDistType onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$mirDistType_List = array('ppm', '100%');" << "\n" << std::flush;
            output << "				$mirDistType_Size = Count( $mirDistType_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $mirDistType_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$mirDistType_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $mirDistType == $mirDistType_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>' . $mirDistType_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Max NU ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $mirDistType == 'ppm' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<input type=text name=MaxHight size=8 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $MaxHight=='' )" << "\n" << std::flush;
            output << "						echo 'MaxHight';" << "\n" << std::flush;
            output << "					else" << "\n" << std::flush;
            output << "						echo $MaxHight;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Multi TSV File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File[] multiple=mutiple hight=>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					foreach( $TSV_File as $TSVs )" << "\n" << std::flush;
            output << "						if( $TSVs == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "							echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='submit' value='Submit' /> " << "\n" << std::flush;
            output << "					</form><br/>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== MirDist ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				//========== svg view var set ==========" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$width = $Top_miRNA * 180;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"<script>" << "\n" << std::flush;
            output << "						var margin = {top: 20, right: 60, bottom: 35, left: 40}," << "\n" << std::flush;
            output << "							width = $width - margin.left - margin.right, //180(per one miRNA)" << "\n" << std::flush;
            output << "							height = 300 - margin.top - margin.bottom;" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var x0 = d3.scale.ordinal()" << "\n" << std::flush;
            output << "							.rangeRoundBands([0, width], .3);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var x1 = d3.scale.ordinal();" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var y = d3.scale.linear()" << "\n" << std::flush;
            output << "							.range([height, 0]);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var color = d3.scale.ordinal()" << "\n" << std::flush;
            output << "							.range(['#F75C2F', '#E8B647', '#838A2D', '#66BAB7', '#6E75A4', '#72636E']);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var xAxis = d3.svg.axis()" << "\n" << std::flush;
            output << "							.scale(x0)" << "\n" << std::flush;
            output << "							.orient('bottom');" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "						var yAxis = d3.svg.axis()" << "\n" << std::flush;
            output << "							.scale(y)" << "\n" << std::flush;
            output << "							.orient('left')" << "\n" << std::flush;
            output << "							.tickFormat(d3.format('.2s'));\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				//==================== svg ====================" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $mirDistType == 'ppm' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					if( $MaxHight == '' || $MaxHight == 'MaxHight' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$MaxValue = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						For( $i = 0; $i < Count( $TSV_File ); ++$i )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$File = $Chart_Type.$TSV_File[$i];" << "\n" << std::flush;
            output << "							$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "							$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								For( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									if( $inFile_Line[$k] >= $MaxValue )" << "\n" << std::flush;
            output << "										$MaxValue = $inFile_Line[$k];" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "					else" << "\n" << std::flush;
            output << "						$MaxValue = $MaxHight;" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					$MaxValue = 100;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < Count( $TSV_File ); ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					$File = $Chart_Type.$TSV_File[$i];" << "\n" << std::flush;
            output << "					$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "					$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$Temp = Tempnam( '/tmp', $TSV_File[$i] );" << "\n" << std::flush;
            output << "					$Ftemp = fopen( $Temp, 'w' );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $mirDistType == 'ppm' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							if( $j >= $Top_miRNA+1 )" << "\n" << std::flush;
            output << "								break;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "					else" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							if( $j >= $Top_miRNA+1 )" << "\n" << std::flush;
            output << "								break;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $j == 0 )" << "\n" << std::flush;
            output << "								fwrite( $Ftemp, $inFile_Lines[$j].\"\\n\" );" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$inFile_Line_clm = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n" << std::flush;
            output << "								$Total_Value = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )" << "\n" << std::flush;
            output << "									$Total_Value += $inFile_Line_clm[ $k ];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								fwrite( $Ftemp, $inFile_Line_clm[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								For( $k = 1; $k < Count( $inFile_Line_clm ); ++$k )" << "\n" << std::flush;
            output << "									fwrite( $Ftemp, \"\\t\".Number_format( $inFile_Line_clm[ $k ]*100/$Total_Value, 0 ));" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								fwrite( $Ftemp, \"\\n\");" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					fclose( $Ftemp );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \"var svg$i = d3.select('body').append('svg')" << "\n" << std::flush;
            output << "							.attr('id', 'svg$i')" << "\n" << std::flush;
            output << "   	 					.attr('width', width + margin.left + margin.right)" << "\n" << std::flush;
            output << "   	 					.attr('height', height + margin.top + margin.bottom)" << "\n" << std::flush;
            output << "   	 					.append('g')" << "\n" << std::flush;
            output << "   	 					.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "						d3.tsv('$Temp', function(error, data) {" << "\n" << std::flush;
            output << "   	 					var sample = d3.keys(data[0]).filter(function(key) { return key !== 'miRNA'; });" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "   	 					data.forEach(function(d) {" << "\n" << std::flush;
            output << "   	 						 d.values = sample.map(function(name) { return {name: name, value: +d[name]}; });" << "\n" << std::flush;
            output << "   	 					});" << "\n" << std::flush;
            output << "   	 		" << "\n" << std::flush;
            output << "   	 					x0.domain(data.map(function(d) { return d.miRNA; }));" << "\n" << std::flush;
            output << "   	 					x1.domain(sample).rangeRoundBands([0, x0.rangeBand()]);" << "\n" << std::flush;
            output << "						y.domain([0, $MaxValue]);" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "							svg$i.append('g')" << "\n" << std::flush;
            output << "   	 						.attr('class', 'x axis')" << "\n" << std::flush;
            output << "   	 						.attr('transform', 'translate(0,' + height + ')')" << "\n" << std::flush;
            output << "   	 						.call(xAxis);" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "							svg$i.append('g')" << "\n" << std::flush;
            output << "   	 						.attr('class', 'y axis')" << "\n" << std::flush;
            output << "   	 						.call(yAxis)" << "\n" << std::flush;
            output << "   	 						.append('text')" << "\n" << std::flush;
            output << "   	 						.attr('transform', 'rotate(-90)')" << "\n" << std::flush;
            output << "   	 						.attr('y', 6)" << "\n" << std::flush;
            output << "   	 						.attr('dy', '.71em')" << "\n" << std::flush;
            output << "   	 						.style('text-anchor', 'end')" << "\n" << std::flush;
            output << "								.text('$TSV_File[$i]');" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "							var barchat = svg$i.selectAll('.barchat')" << "\n" << std::flush;
            output << "   	 						.data(data)" << "\n" << std::flush;
            output << "   	 						.enter().append('g')" << "\n" << std::flush;
            output << "   	 						.attr('class', 'g')" << "\n" << std::flush;
            output << "   	 						.attr('transform', function(d) { return 'translate(' + x0(d.miRNA) + ',0)'; });" << "\n" << std::flush;
            output << "   	 		" << "\n" << std::flush;
            output << "   	 					var div = d3.select('body').append('div')" << "\n" << std::flush;
            output << "   	 						.attr('class', 'tooltip')" << "\n" << std::flush;
            output << "							.attr('id', 'Mir')" << "\n" << std::flush;
            output << "   	 						.style('opacity', 0);" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "   	 					barchat.selectAll('rect')" << "\n" << std::flush;
            output << "   	 						.data(function(d) { return d.values; })" << "\n" << std::flush;
            output << "   	 						.enter().append('rect')" << "\n" << std::flush;
            output << "   	 						.attr('width', x1.rangeBand() - 5)" << "\n" << std::flush;
            output << "   	 						.attr('x', function(d) { return x1(d.name); })" << "\n" << std::flush;
            output << "   	 						.attr('y', function(d){" << "\n" << std::flush;
            output << "   	 							if( d.value < $MaxValue )" << "\n" << std::flush;
            output << "   	 								return y(d.value);" << "\n" << std::flush;
            output << "   	 							else" << "\n" << std::flush;
            output << "   	 								return y($MaxValue);" << "\n" << std::flush;
            output << "   	 						})" << "\n" << std::flush;
            output << "   	 						.attr('height', function(d){" << "\n" << std::flush;
            output << "   	 							if( d.value < $MaxValue )" << "\n" << std::flush;
            output << "   	 								return height - y(d.value);" << "\n" << std::flush;
            output << "   	 							else" << "\n" << std::flush;
            output << "   	 								return height - y($MaxValue);" << "\n" << std::flush;
            output << "   	 						})" << "\n" << std::flush;
            output << "   	 						.style('fill', function(d) { return color(d.name); })" << "\n" << std::flush;
            output << "   	 						.on('mouseover', function(d) {" << "\n" << std::flush;
            output << "   	 							div.transition()" << "\n" << std::flush;
            output << "   	 								.duration(200)" << "\n" << std::flush;
            output << "   	 								.style('opacity', .9);" << "\n" << std::flush;
            output << "   	 		" << "\n" << std::flush;
            output << "								div.html('<table align=center ><tr><th>Sample</th><th>Value</th></tr><tr><th>' +" << "\n" << std::flush;
            output << "										 d.name + '</th><th>' + d.value + '</th></tr>')" << "\n" << std::flush;
            output << "   	 								.style('left', (d3.event.pageX) + 'px')" << "\n" << std::flush;
            output << "   	 								.style('top', (d3.event.pageY - 28) + 'px');" << "\n" << std::flush;
            output << "   	 						})" << "\n" << std::flush;
            output << "   	 						.on('mouseout', function(d) {" << "\n" << std::flush;
            output << "   	 							div.transition()" << "\n" << std::flush;
            output << "   	 								.duration(500)" << "\n" << std::flush;
            output << "   	 								.style('opacity', 0);" << "\n" << std::flush;
            output << "   	 						});" << "\n" << std::flush;
            output << "   	 		" << "\n" << std::flush;
            output << "							var legend = svg$i.selectAll('.legend')" << "\n" << std::flush;
            output << "   	 						.data(sample.slice())" << "\n" << std::flush;
            output << "   	 						.enter().append('g')" << "\n" << std::flush;
            output << "   	 						.attr('class', 'legend')" << "\n" << std::flush;
            output << "   	 						.attr('transform', function(d, i) { return 'translate(0,' + i * 20 + ')'; });" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "   	 					legend.append('rect')" << "\n" << std::flush;
            output << "   	 						.attr('x', width + 40)" << "\n" << std::flush;
            output << "   	 						.attr('width', 18)" << "\n" << std::flush;
            output << "   	 						.attr('height', 18)" << "\n" << std::flush;
            output << "   	 						.style('fill', color);" << "\n" << std::flush;
            output << "   	 				" << "\n" << std::flush;
            output << "   	 					legend.append('text')" << "\n" << std::flush;
            output << "   	 						.attr('x', width + 34)" << "\n" << std::flush;
            output << "   	 						.attr('y', 9)" << "\n" << std::flush;
            output << "   	 						.attr('dy', '.35em')" << "\n" << std::flush;
            output << "   	 						.style('text-anchor', 'end')" << "\n" << std::flush;
            output << "   	 						.text(function(d) { return d; });" << "\n" << std::flush;
            output << "   	 				});\";" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo '</script>';" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type == 'DotPlot/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== TSV File ====================-->" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $TSV_File == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo '</select>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File2 onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($TSV_File2=='') echo 'selected'; echo '>Select TSV</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $TSV_File2 == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== is_Abundant ====================-->" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				$isAbundant_List = array('MostAbundant', 'AllmiRNA');" << "\n" << std::flush;
            output << "				$isAbundant_Size = Count( $isAbundant_List );" << "\n" << std::flush;
            output << "				if( $isAbundant == '' )" << "\n" << std::flush;
            output << "					$isAbundant = 'MostAbundant';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$isAbundant_List[$i].' ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					if( $isAbundant == $isAbundant_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					echo '>' . $isAbundant_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== isLog ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=isLog onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($isLog=='') echo 'selected'; echo 'value= >isLog</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$isLog_List = array(2, 4, 6, 8, 10);" << "\n" << std::flush;
            output << "				$isLog_Size = Count( $isLog_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $isLog_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$isLog_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $isLog == $isLog_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$isLog_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMin' value='$ForceMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceMax' value='$ForceMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== ForceMin&Max ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<input type=text name=ForceMin size=5 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceMin=='' )" << "\n" << std::flush;
            output << "					echo 'ForceMin';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $ForceMin;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "				echo '<input type=text name=ForceMax size=5 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceMax=='' )" << "\n" << std::flush;
            output << "					echo 'ForceMax';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $ForceMax;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter NU ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<input type=text name=Filter size=4 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $Filter=='' )" << "\n" << std::flush;
            output << "					echo 'Filter';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $Filter;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter GMPM ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=FGMPM >';" << "\n" << std::flush;
            output << "				echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$FGMPM_List = array('GMPM', 'GM', 'PM', 'Tailing');" << "\n" << std::flush;
            output << "				$FGMPM_Size = Count( $FGMPM_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $FGMPM_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$FGMPM_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $FGMPM == $FGMPM_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>' . $FGMPM_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='GMPM' value='$GMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Bar_Pie' value='$Bar_Pie' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File2' value='$TSV_File2' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "					<input type='submit' value='Submit' /> " << "\n" << std::flush;
            output << "					</form><br/>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Index & Header ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$File1 = $Chart_Type.$TSV_File; 	//y" << "\n" << std::flush;
            output << "				$File2 = $Chart_Type.$TSV_File2;	//x" << "\n" << std::flush;
            output << "				$Files = array( $File1, $File2 );" << "\n" << std::flush;
            output << "				$Index = array();" << "\n" << std::flush;
            output << "				$Column = 0;" << "\n" << std::flush;
            output << "				$FColumn = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $isAbundant != 'MostAbundant' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					For( $l = 0; $l < Count( $Files ); ++$l )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$inFile = File_get_contents( $Files[$l] );" << "\n" << std::flush;
            output << "						$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $i == 0 )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								For( $j = 0; $j < Count( $inFile_Line ); ++$j )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									if( $inFile_Line[$j] == $GMPM )" << "\n" << std::flush;
            output << "										$Column = $j;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "									if( $inFile_Line[$j] == $FGMPM )" << "\n" << std::flush;
            output << "										$FColumn = $j;" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$miRNA_Seed = Explode( '*', $inFile_Line[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( Count( $miRNA_Seed ) != 2 )" << "\n" << std::flush;
            output << "									Array_Push( $Index, $inFile_Line[0] );" << "\n" << std::flush;
            output << "								else" << "\n" << std::flush;
            output << "									Array_Push( $Index, $miRNA_Seed[1] );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					For( $l = 0; $l < Count( $Files ); ++$l )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$inFile = File_get_contents( $Files[$l] );" << "\n" << std::flush;
            output << "						$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $i == 0 )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								For( $j = 0; $j < Count( $inFile_Line ); ++$j )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									if( $inFile_Line[$j] == $GMPM )" << "\n" << std::flush;
            output << "										$Column = $j;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "									if( $inFile_Line[$j] == $FGMPM )" << "\n" << std::flush;
            output << "										$FColumn = $j;" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								Array_Push( $Index, 'miRNA' );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$miRNA_Seed = Explode( '*', $inFile_Line[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( Count( $miRNA_Seed ) != 2 )" << "\n" << std::flush;
            output << "								{}" << "\n" << std::flush;
            output << "								else" << "\n" << std::flush;
            output << "									Array_Push( $Index, $miRNA_Seed[1] );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Sort( $Index, SORT_STRING );" << "\n" << std::flush;
            output << "				$uIndex = Array_Unique( $Index, SORT_STRING );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Read File & Log ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$Anno_Value1 = array();" << "\n" << std::flush;
            output << "				$Anno_Value2 = array();" << "\n" << std::flush;
            output << "				$Filter_Array1 = array();" << "\n" << std::flush;
            output << "				$Filter_Array2 = array();" << "\n" << std::flush;
            output << "				$Sample_Name = array();" << "\n" << std::flush;
            output << "				$FSample_Name = array();" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $l = 0; $l < Count( $Files ); ++$l )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					$inFile = File_get_contents( $Files[$l] );" << "\n" << std::flush;
            output << "					$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $i == 0 )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							if( $l == 0 )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$Anno_Value1['miRNA'] = $inFile_Line[0];" << "\n" << std::flush;
            output << "								$Filter_Array1['miRNA']=$inFile_Line[0].'_F'.$FGMPM;" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$Anno_Value2['miRNA'] = $inFile_Line[0];" << "\n" << std::flush;
            output << "								$Filter_Array2['miRNA']=$inFile_Line[0].'_F'.$FGMPM;" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							Array_Push( $Sample_Name, $inFile_Line[0] );" << "\n" << std::flush;
            output << "							Array_Push( $FSample_Name, $inFile_Line[0].'_F'.$FGMPM );" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "						else" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$miRNA_Seed = Explode( '*', $inFile_Line[0] );" << "\n" << std::flush;
            output << "							$Value = $inFile_Line[$Column];" << "\n" << std::flush;
            output << "							$FValue = $inFile_Line[$FColumn];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $isLog != '' && $Value != 0 )" << "\n" << std::flush;
            output << "								$Value = ( Log($inFile_Line[$Column]) / Log($isLog) );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $isLog != '' && $FValue != 0 )" << "\n" << std::flush;
            output << "								$FValue = ( Log($inFile_Line[$FColumn]) / Log($isLog) );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( Count( $miRNA_Seed ) != 2 )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								if( $l == 0 )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									$Anno_Value1[$inFile_Line[0]] = Round($Value,2);" << "\n" << std::flush;
            output << "									$Filter_Array1[$inFile_Line[0]]=Round($FValue,2);" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "								else" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									$Anno_Value2[$inFile_Line[0]] = Round($Value,2);" << "\n" << std::flush;
            output << "									$Filter_Array2[$inFile_Line[0]]=Round($FValue,2);" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								if( $l == 0 )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									$Anno_Value1[$miRNA_Seed[1]] = Round($Value,2);" << "\n" << std::flush;
            output << "									$Filter_Array1[$miRNA_Seed[1]]=Round($FValue,2);" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "								else" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									$Anno_Value2[$miRNA_Seed[1]] = Round($Value,2);" << "\n" << std::flush;
            output << "									$Filter_Array2[$miRNA_Seed[1]]=Round($FValue,2);" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter & Temp ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$Temp = Tempnam( '/tmp', $GMPM.'_'.$Sample_Name[0].'_'.$Sample_Name[1].'_'.$isAbundant.'_'.$isLog.'_'.$Filter.'_'.$FGMPM );" << "\n" << std::flush;
            output << "				$Ftemp = Fopen( $Temp, 'w' );" << "\n" << std::flush;
            output << "				$MaxAxis = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Fwrite( $Ftemp, 'miRNA'.\"\\t\"." << "\n" << std::flush;
            output << "						$Anno_Value1['miRNA'].\"\\t\"." << "\n" << std::flush;
            output << "						$Anno_Value2['miRNA'].\"\\t\"." << "\n" << std::flush;
            output << "						$Filter_Array1['miRNA'].\"\\t\"." << "\n" << std::flush;
            output << "						$Filter_Array2['miRNA'].\"\\n\" );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < Count( $Index ); ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					if( $uIndex[$i] != '' && $uIndex[$i] != 'miRNA' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						if( $Anno_Value1[$uIndex[$i]] == '' )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$Anno_Value1[$uIndex[$i]]  = 0;" << "\n" << std::flush;
            output << "							$Filter_Array1[$uIndex[$i]]= 0;" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $Anno_Value2[$uIndex[$i]] == '' )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$Anno_Value2[$uIndex[$i]]  = 0;" << "\n" << std::flush;
            output << "							$Filter_Array2[$uIndex[$i]]= 0;" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $Filter != 'Filter' && $FGMPM != '' )" << "\n" << std::flush;
            output << "							if( $Filter_Array1[$uIndex[$i]] < $Filter && $Filter_Array2[$uIndex[$i]] < $Filter )" << "\n" << std::flush;
            output << "								continue;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						Fwrite( $Ftemp, $uIndex[$i].\"\\t\"." << "\n" << std::flush;
            output << "										$Anno_Value1[$uIndex[$i]].\"\\t\"." << "\n" << std::flush;
            output << "										$Anno_Value2[$uIndex[$i]].\"\\t\"." << "\n" << std::flush;
            output << "										$Filter_Array1[$uIndex[$i]].\"\\t\"." << "\n" << std::flush;
            output << "										$Filter_Array2[$uIndex[$i]].\"\\n\" );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $MaxAxis < $Anno_Value1[$uIndex[$i]] )" << "\n" << std::flush;
            output << "							$MaxAxis = $Anno_Value1[$uIndex[$i]];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $MaxAxis < $Anno_Value2[$uIndex[$i]] )" << "\n" << std::flush;
            output << "							$MaxAxis = $Anno_Value2[$uIndex[$i]];" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Fclose( $Ftemp );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $FGMPM == '' )" << "\n" << std::flush;
            output << "					$FGMPM = 'Filter';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceMin == 'ForceMin' || $ForceMin == '' )" << "\n" << std::flush;
            output << "					$ForceMin = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceMax == 'ForceMax' || $ForceMax == '' )" << "\n" << std::flush;
            output << "					$ForceMax = $MaxAxis;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== DotPlot ====================-->" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "				echo \"<script>" << "\n" << std::flush;
            output << "					var svg_width  = window.innerWidth;" << "\n" << std::flush;
            output << "					var svg_height = window.innerHeight;" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var margin = {top: 20, right: 60, bottom: 60, left: 60}" << "\n" << std::flush;
            output << "						width = svg_width - margin.left - margin.right," << "\n" << std::flush;
            output << "						height = svg_height - margin.top - margin.bottom;" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var x = d3.scale.linear()" << "\n" << std::flush;
            output << "						.range([0, width]);" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var y = d3.scale.linear()" << "\n" << std::flush;
            output << "						.range([height, 0]);" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var color = d3.scale.category10();" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var xAxis = d3.svg.axis()" << "\n" << std::flush;
            output << "						.scale(x)" << "\n" << std::flush;
            output << "						.orient('bottom');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var yAxis = d3.svg.axis()" << "\n" << std::flush;
            output << "						.scale(y)" << "\n" << std::flush;
            output << "						.orient('left');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					var svg = d3.select('body').append('svg')" << "\n" << std::flush;
            output << "						.attr('width', width)" << "\n" << std::flush;
            output << "						.attr('height', height)" << "\n" << std::flush;
            output << "					.append('g')" << "\n" << std::flush;
            output << "						.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "					d3.tsv('$Temp', function(error, data) {" << "\n" << std::flush;
            output << "						data.forEach(function(d) {" << "\n" << std::flush;
            output << "							d.$Sample_Name[1] = +d.$Sample_Name[1];" << "\n" << std::flush;
            output << "							d.$Sample_Name[0] = +d.$Sample_Name[0];" << "\n" << std::flush;
            output << "						});" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						var xExtent = d3.extent(data, function(d) { return d.$Sample_Name[1]; });" << "\n" << std::flush;
            output << "						var yExtent = d3.extent(data, function(d) { return d.$Sample_Name[0]; });" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						xExtent[0] = $ForceMin;" << "\n" << std::flush;
            output << "						yExtent[0] = $ForceMin;" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						xExtent[1] = $ForceMax;" << "\n" << std::flush;
            output << "						yExtent[1] = $ForceMax;" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						x.domain( xExtent ).nice();" << "\n" << std::flush;
            output << "						y.domain( yExtent ).nice();" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						svg.append('g')" << "\n" << std::flush;
            output << "							.attr('class', 'x axis')" << "\n" << std::flush;
            output << "							.attr('transform', 'translate(0,' + height + ')')" << "\n" << std::flush;
            output << "							.call(xAxis)" << "\n" << std::flush;
            output << "							.append('text')" << "\n" << std::flush;
            output << "							.attr('class', 'label')" << "\n" << std::flush;
            output << "							.attr('x', width)" << "\n" << std::flush;
            output << "							.attr('y', -6)" << "\n" << std::flush;
            output << "							.style('text-anchor', 'end')" << "\n" << std::flush;
            output << "							.text('$Sample_Name[1]');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						svg.append('g')" << "\n" << std::flush;
            output << "							.attr('class', 'y axis')" << "\n" << std::flush;
            output << "							.call(yAxis)" << "\n" << std::flush;
            output << "							.append('text')" << "\n" << std::flush;
            output << "							.attr('class', 'label')" << "\n" << std::flush;
            output << "							.attr('transform', 'rotate(-90)')" << "\n" << std::flush;
            output << "							.attr('y', 6)" << "\n" << std::flush;
            output << "							.attr('dy', '.71em')" << "\n" << std::flush;
            output << "							.style('text-anchor', 'end')" << "\n" << std::flush;
            output << "							.text('$Sample_Name[0]');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						svg.append('g')" << "\n" << std::flush;
            output << "							.attr('class', 'xyline')" << "\n" << std::flush;
            output << "							.append('line')" << "\n" << std::flush;
            output << "							.attr('x1', x($ForceMin) )" << "\n" << std::flush;
            output << "							.attr('y1', y($ForceMin) )" << "\n" << std::flush;
            output << "							.attr('x2', x($ForceMax) )" << "\n" << std::flush;
            output << "							.attr('y2', y($ForceMax) )" << "\n" << std::flush;
            output << "							.style('stroke','rgb(255,0,0)')" << "\n" << std::flush;
            output << "							.style('stroke-width','1');" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						var div = d3.select('body').append('div')" << "\n" << std::flush;
            output << "							.attr('class', 'tooltip')" << "\n" << std::flush;
            output << "							.attr('id', 'Dot')" << "\n" << std::flush;
            output << "							.style('opacity', 0);" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						svg.selectAll('.dot')" << "\n" << std::flush;
            output << "							.data(data)" << "\n" << std::flush;
            output << "							.enter().append('circle')" << "\n" << std::flush;
            output << "							.attr('class', 'dot')" << "\n" << std::flush;
            output << "							.attr('r', 5)" << "\n" << std::flush;
            output << "							.attr('cx', function(d) { return x(d.$Sample_Name[1]); })" << "\n" << std::flush;
            output << "							.attr('cy', function(d) { return y(d.$Sample_Name[0]); })" << "\n" << std::flush;
            output << "							.style('fill', function(d) { return color(d.species); })" << "\n" << std::flush;
            output << "							.on('mouseover', function(d) {" << "\n" << std::flush;
            output << "								div.transition()" << "\n" << std::flush;
            output << "									.duration(200)" << "\n" << std::flush;
            output << "									.style('opacity', .9);" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "								div.html('<table align=center ><tr><th colspan=3 >' + d.miRNA +" << "\n" << std::flush;
            output << "										 '</th></tr><tr><th>Sample</th><th>$GMPM</th><th>$FGMPM</th></tr><tr><th>$Sample_Name[0]</th><th>' +" << "\n" << std::flush;
            output << "										 d.$Sample_Name[0] + '</th><th>' + d.$Sample_Name[0]_F$FGMPM + '</th></tr><tr><th>$Sample_Name[1]</th><th>' +" << "\n" << std::flush;
            output << "										 d.$Sample_Name[1] + '</th><th>' + d.$Sample_Name[1]_F$FGMPM + '</th></tr>')" << "\n" << std::flush;
            output << "									.style('left', (d3.event.pageX) + 'px')" << "\n" << std::flush;
            output << "									.style('top', (d3.event.pageY - 28) + 'px');" << "\n" << std::flush;
            output << "							})" << "\n" << std::flush;
            output << "							.on('mouseout', function(d) {" << "\n" << std::flush;
            output << "								div.transition()" << "\n" << std::flush;
            output << "									.duration(500)" << "\n" << std::flush;
            output << "									.style('opacity', 0);" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "							window.onresize = function(){" << "\n" << std::flush;
            output << "							};" << "\n" << std::flush;
            output << "					});" << "\n" << std::flush;
            output << "				</script>\";" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "		}" << "\n" << std::flush;
            output << "		else" << "\n" << std::flush;
            output << "		{" << "\n" << std::flush;
            output << "			if( $Chart_Type == 'ValPlot/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '	<script src=http://140.113.239.177/joye/ForAgoSorting/lib/d3.min.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/jquery-1.7.min.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/jquery.event.drag-2.0.min.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/slick.core.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/slick.grid.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/slick.dataview.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/d3.parcoords.js ></script>" << "\n" << std::flush;
            output << "						<script src=http://140.113.239.177/joye/ForAgoSorting/lib/divgrid.js ></script>" << "\n" << std::flush;
            output << "						<link rel=stylesheet type=text/css href=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/slick.grid.css />" << "\n" << std::flush;
            output << "						<link rel=stylesheet type=text/css href=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/jquery-ui-1.8.16.custom.css />" << "\n" << std::flush;
            output << "						<link rel=stylesheet type=text/css href=http://140.113.239.177/joye/ForAgoSorting/lib/slickgrid/examples.css />" << "\n" << std::flush;
            output << "						<link rel=stylesheet type=text/css href=http://140.113.239.177/joye/ForAgoSorting/lib/d3.parcoords.css />" << "\n" << std::flush;
            output << "						<link rel=stylesheet type=text/css href=http://140.113.239.177/joye/ForAgoSorting/lib/style.css />';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== TSV File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $TSV_File == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== is_Abundant ====================-->" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				$isAbundant_List = array('MostAbundant', 'AllmiRNA');" << "\n" << std::flush;
            output << "				$isAbundant_Size = Count( $isAbundant_List );" << "\n" << std::flush;
            output << "				if( $isAbundant == '' )" << "\n" << std::flush;
            output << "					$isAbundant = 'MostAbundant';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$isAbundant_List[$i].' ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					if( $isAbundant == $isAbundant_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					echo '>' . $isAbundant_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== isLog ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=isLog onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($isLog=='') echo 'selected'; echo 'value= >isLog</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$isLog_List = array(2, 4, 6, 8, 10);" << "\n" << std::flush;
            output << "				$isLog_Size = Count( $isLog_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $isLog_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$isLog_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $isLog == $isLog_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$isLog_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='FGMPM' value='$FGMPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Filter' value='$Filter' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMin' value='$FilterMin' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='FilterMax' value='$FilterMax' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter Min & Max ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<input type=text name=FilterMin size=5 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $FilterMin=='' )" << "\n" << std::flush;
            output << "					echo 'FilterMin';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $FilterMin;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "				echo '<input type=text name=FilterMax size=5 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $FilterMax=='' )" << "\n" << std::flush;
            output << "					echo 'FilterMax';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $FilterMax;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter NU ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<input type=text name=Filter size=4 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $Filter=='' )" << "\n" << std::flush;
            output << "					echo 'Filter';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $Filter;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Filter GMPM ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=FGMPM >';" << "\n" << std::flush;
            output << "				echo '<option '; if($FGMPM=='') echo 'selected'; echo ' value= >GM or PM</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $FGMPM == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"" << "\n" << std::flush;
            output << "					<input type='hidden' name='isLog' value='$isLog' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='isAbundant' value='$isAbundant' />" << "\n" << std::flush;
            output << "					<input type='submit' value='Submit' /> " << "\n" << std::flush;
            output << "					</form><br/>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Read File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $TSV_File != '' )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					if( $FGMPM != '' && $Filter != '' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$FFile = $Chart_Type.$FGMPM;" << "\n" << std::flush;
            output << "						$FinFile = File_get_contents( $FFile );" << "\n" << std::flush;
            output << "						$FinFile_Lines = Explode( \"\\n\", $FinFile );" << "\n" << std::flush;
            output << "						$FArray = Array();" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$File = $Chart_Type.$TSV_File;" << "\n" << std::flush;
            output << "					$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "					$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$Temp = Tempnam( '/tmp', $TSV_File );" << "\n" << std::flush;
            output << "					$Ftemp = Fopen( $Temp, 'w' );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					$Sample_Nu = Count( Explode( \"\\t\", $inFile_Lines[0] ))-1;" << "\n" << std::flush;
            output << "					$Value_Max = 0;" << "\n" << std::flush;
            output << "					$Value_Min = 0;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					For( $i = 0; $i < Count( $inFile_Lines )-1; ++$i )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						if( $FGMPM != '' && $Filter != '' && $i != 0 )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$Check_Filter = Array();" << "\n" << std::flush;
            output << "							$FinFile_Line = Explode( \"\\t\", $FinFile_Lines[$i] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							For( $j = 1; $j < Count( $FinFile_Line ); ++$j )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$FValue = $FinFile_Line[$j];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( $isLog != '' && $FValue != 0 )" << "\n" << std::flush;
            output << "									$FValue = ( Log($FinFile_Line[$j]) / Log($isLog) );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( $FValue <= $Filter )" << "\n" << std::flush;
            output << "									Array_Push( $Check_Filter, $FValue );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( Count( $Check_Filter ) == Count( $FinFile_Line )-1 )" << "\n" << std::flush;
            output << "								Continue;" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						$Check_Filter = 'keep';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $i != 0 )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							$Check_Min = Array();" << "\n" << std::flush;
            output << "							$Check_Max = Array();" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							$inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$Value = $inFile_Line[$j];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( $isLog != '' && $Value != 0 )" << "\n" << std::flush;
            output << "									$Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( $FilterMin != 'FilterMin' && $FilterMin != '' )" << "\n" << std::flush;
            output << "									if( $Value <= $FilterMin )" << "\n" << std::flush;
            output << "										Array_Push( $Check_Min, $Value );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( $FilterMax != 'FilterMax' && $FilterMax != '' )" << "\n" << std::flush;
            output << "									if( $Value >= $FilterMax )" << "\n" << std::flush;
            output << "										Array_Push( $Check_Max, $Value );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( Count( $Check_Min ) == Count( $inFile_Line )-1 || Count( $Check_Max ) == Count( $inFile_Line )-1 )" << "\n" << std::flush;
            output << "								$Check_Filter = 'notkeep';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $isAbundant == 'MostAbundant' )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$miRNA = Explode( '*', $inFile_Line[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								if( Count( $miRNA ) != 2 )" << "\n" << std::flush;
            output << "									$Check_Filter = 'notkeep';" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "							if( $Check_Filter == 'keep' )" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									$Value = $inFile_Line[$j];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "									if( $isLog != '' && $Value != 0 )" << "\n" << std::flush;
            output << "										$Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "									if( $Value_Max <= $Value )" << "\n" << std::flush;
            output << "										$Value_Max =  $Value;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "									if( $Value_Min >= $Value )" << "\n" << std::flush;
            output << "										$Value_Min =  $Value;" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $Check_Filter == 'keep' )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							if( $isLog != '' && $i != 0)" << "\n" << std::flush;
            output << "							{" << "\n" << std::flush;
            output << "								$inFile_Line = Explode( \"\\t\", $inFile_Lines[$i] );" << "\n" << std::flush;
            output << "								Fwrite( $Ftemp, $inFile_Line[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								For( $j = 1; $j < Count( $inFile_Line ); ++$j )" << "\n" << std::flush;
            output << "								{" << "\n" << std::flush;
            output << "									if( $inFile_Line[$j] != 0 )" << "\n" << std::flush;
            output << "									{	" << "\n" << std::flush;
            output << "										$Value = ( Log($inFile_Line[$j]) / Log($isLog) );" << "\n" << std::flush;
            output << "										Fwrite( $Ftemp, \"\\t\".$Value );" << "\n" << std::flush;
            output << "									}" << "\n" << std::flush;
            output << "									else" << "\n" << std::flush;
            output << "										Fwrite( $Ftemp, \"\\t\".'0' );" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "								Fwrite( $Ftemp, \"\\n\" );" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "							else" << "\n" << std::flush;
            output << "								Fwrite( $Ftemp, $inFile_Lines[$i].\"\\n\" );" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, 'MinValue' );" << "\n" << std::flush;
            output << "					For( $i = 0; $i < $Sample_Nu; ++$i )" << "\n" << std::flush;
            output << "						Fwrite( $Ftemp, \"\\t\".$Value_Min );" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, \"\\n\" );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, 'MaxValue' );" << "\n" << std::flush;
            output << "					For( $i = 0; $i < $Sample_Nu; ++$i )" << "\n" << std::flush;
            output << "						Fwrite( $Ftemp, \"\\t\".$Value_Max );" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, \"\\n\" );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					fclose( $Ftemp );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== ValPlot ====================-->" << "\n" << std::flush;
            output << "	" << "\n" << std::flush;
            output << "					echo \"<div id='example' class='parcoords' style='height:240px;'></div>" << "\n" << std::flush;
            output << "						<div id='grid'></div>" << "\n" << std::flush;
            output << "						<script id='brushing'>" << "\n" << std::flush;
            output << "							var parcoords = d3.parcoords()('#example')" << "\n" << std::flush;
            output << "									.alpha(0.4)" << "\n" << std::flush;
            output << "									.mode('queue') // progressive rendering" << "\n" << std::flush;
            output << "									.height(d3.max([document.body.clientHeight-326, 220]))" << "\n" << std::flush;
            output << "									.margin({" << "\n" << std::flush;
            output << "										top: 36," << "\n" << std::flush;
            output << "										left: 0," << "\n" << std::flush;
            output << "										right: 0," << "\n" << std::flush;
            output << "										bottom: 16" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "							// load csv file and create the chart" << "\n" << std::flush;
            output << "							d3.tsv('$Temp', function(data) {" << "\n" << std::flush;
            output << "								// slickgrid needs each data element to have an id" << "\n" << std::flush;
            output << "								data.forEach(function(d,i) { d.id = d.id || i; });" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								parcoords" << "\n" << std::flush;
            output << "									.data(data)" << "\n" << std::flush;
            output << "									.hideAxis(['miRNA'])" << "\n" << std::flush;
            output << "									.hideAxis(['id'])" << "\n" << std::flush;
            output << "									.render()" << "\n" << std::flush;
            output << "									.reorderable()" << "\n" << std::flush;
            output << "									.brushMode('1D-axes');" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// setting up grid" << "\n" << std::flush;
            output << "								var column_keys = d3.keys(data[0]);" << "\n" << std::flush;
            output << "								var columns = column_keys.map(function(key,i) {" << "\n" << std::flush;
            output << "									return {" << "\n" << std::flush;
            output << "										id: key," << "\n" << std::flush;
            output << "										name: key," << "\n" << std::flush;
            output << "										field: key," << "\n" << std::flush;
            output << "										sortable: true" << "\n" << std::flush;
            output << "									}" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								var options = {" << "\n" << std::flush;
            output << "									enableCellNavigation: true," << "\n" << std::flush;
            output << "									enableColumnReorder: false," << "\n" << std::flush;
            output << "									multiColumnSort: false" << "\n" << std::flush;
            output << "								};" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								var dataView = new Slick.Data.DataView();" << "\n" << std::flush;
            output << "								var grid = new Slick.Grid('#grid', dataView, columns, options);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// wire up model events to drive the grid" << "\n" << std::flush;
            output << "								dataView.onRowCountChanged.subscribe(function (e, args) {" << "\n" << std::flush;
            output << "									grid.updateRowCount();" << "\n" << std::flush;
            output << "									grid.render();" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								dataView.onRowsChanged.subscribe(function (e, args) {" << "\n" << std::flush;
            output << "									grid.invalidateRows(args.rows);" << "\n" << std::flush;
            output << "									grid.render();" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// column sorting" << "\n" << std::flush;
            output << "								var sortcol = column_keys[0];" << "\n" << std::flush;
            output << "								var sortdir = 1;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								function comparer(a, b) {" << "\n" << std::flush;
            output << "									var x = a[sortcol], y = b[sortcol];" << "\n" << std::flush;
            output << "									return (x == y ? 0 : (x > y ? 1 : -1));" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "								" << "\n" << std::flush;
            output << "								// click header to sort grid column" << "\n" << std::flush;
            output << "								grid.onSort.subscribe(function (e, args) {" << "\n" << std::flush;
            output << "									sortdir = args.sortAsc ? 1 : -1;" << "\n" << std::flush;
            output << "									sortcol = args.sortCol.field;" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "									if ($.browser.msie && $.browser.version <= 8) {" << "\n" << std::flush;
            output << "										dataView.fastSort(sortcol, args.sortAsc);" << "\n" << std::flush;
            output << "									} else {" << "\n" << std::flush;
            output << "										dataView.sort(comparer, args.sortAsc);" << "\n" << std::flush;
            output << "									}" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// highlight row in chart" << "\n" << std::flush;
            output << "								grid.onMouseEnter.subscribe(function(e,args) {" << "\n" << std::flush;
            output << "									var i = grid.getCellFromEvent(e).row;" << "\n" << std::flush;
            output << "									var d = parcoords.brushed() || data;" << "\n" << std::flush;
            output << "									parcoords.highlight([d[i]]);" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "								grid.onMouseLeave.subscribe(function(e,args) {" << "\n" << std::flush;
            output << "									parcoords.unhighlight();" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// fill grid with data" << "\n" << std::flush;
            output << "								gridUpdate(data);" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								// update grid on brush" << "\n" << std::flush;
            output << "								parcoords.on('brush', function(d) {" << "\n" << std::flush;
            output << "									gridUpdate(d);" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "							" << "\n" << std::flush;
            output << "								function gridUpdate(data) {" << "\n" << std::flush;
            output << "									dataView.beginUpdate();" << "\n" << std::flush;
            output << "									dataView.setItems(data);" << "\n" << std::flush;
            output << "									dataView.endUpdate();" << "\n" << std::flush;
            output << "								};" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "						</script>\";" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type == 'LenPlus/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.js></script>';" << "\n" << std::flush;
            output << "				echo '<link href=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';" << "\n" << std::flush;
            output << "				echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/head.min.js ></script>';" << "\n" << std::flush;
            output << "				echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/reveal.min.js></script>';" << "\n" << std::flush;
            output << "				echo '<link rel=stylesheet href=http://140.113.239.177/joye/ForAgoSorting/lib/reveal.min.css>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== ReadCount & PPM ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<select name=RCPPM onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$RCPPM_List = array('ppm', 'readcount');" << "\n" << std::flush;
            output << "				$RCPPM_Size = Count( $RCPPM_List );" << "\n" << std::flush;
            output << "				if( $RCPPM == '' )" << "\n" << std::flush;
            output << "					$RCPPM = 'ppm';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $RCPPM_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$RCPPM_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $RCPPM == $RCPPM_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>' . $RCPPM_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceY' value='$ForceY' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== ForceY ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<input type=text name=ForceY size=3 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceY=='' )" << "\n" << std::flush;
            output << "					echo 'Hight';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $ForceY;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='RCPPM' value='$RCPPM' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='submit' value='Submit' /> " << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== TSV File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep '.$RCPPM.'.tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<div class=reveal><div class=slides>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					$AGO = Explode( '_', $TSV_List[$i]);" << "\n" << std::flush;
            output << "					echo '<section>';" << "\n" << std::flush;
            output << "					echo '<h2>'.$AGO[0].'</h2>';" << "\n" << std::flush;
            output << "					echo '<svg id=bar'.$i.'></svg>';" << "\n" << std::flush;
            output << "					echo '</section>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Slide Set ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '</div></div>';" << "\n" << std::flush;
            output << "				echo \"<script>" << "\n" << std::flush;
            output << "						Reveal.initialize({" << "\n" << std::flush;
            output << "							center: false," << "\n" << std::flush;
            output << "							transition: Reveal.getQueryHash().transition || 'fade', // default/cube/page/concave/zoom/linear/fade/none" << "\n" << std::flush;
            output << "						});" << "\n" << std::flush;
            output << "					</script>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== Len Chart ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					$File = $Chart_Type.$TSV_List[$i];" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \"<script>" << "\n" << std::flush;
            output << "							d3.tsv( '$File', function( tsv_data ) {" << "\n" << std::flush;
            output << "							var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'Length'; });" << "\n" << std::flush;
            output << "							var data = new Array();" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "							tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "								d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "									return {" << "\n" << std::flush;
            output << "										x: Number(key)," << "\n" << std::flush;
            output << "										y: +d[key]" << "\n" << std::flush;
            output << "									};" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "								});" << "\n" << std::flush;
            output << "								var sample = new Object();" << "\n" << std::flush;
            output << "								sample.key = d['Length'];" << "\n" << std::flush;
            output << "								sample.values = d.len;" << "\n" << std::flush;
            output << "								data.push( sample );" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "							nv.addGraph({" << "\n" << std::flush;
            output << "								generate: function() {" << "\n" << std::flush;
            output << "									var width = nv.utils.windowSize().width," << "\n" << std::flush;
            output << "										height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "									var chart = nv.models.multiBarChart()\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $ForceY != '' && $ForceY != 'Hight' )" << "\n" << std::flush;
            output << "							echo \"		.forceY([$ForceY,0])\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo \"				.width(width)" << "\n" << std::flush;
            output << "										.height(height)" << "\n" << std::flush;
            output << "										.color( ['#000000', '#FF0000', '#0000FF', '#FFBF00', '#088A08', '#6E6E6E'])" << "\n" << std::flush;
            output << "										.stacked(true);" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "									var svg = d3.select('#bar$i').datum(data);" << "\n" << std::flush;
            output << "									console.log('calling chart');" << "\n" << std::flush;
            output << "									svg.transition().duration(0).call(chart);" << "\n" << std::flush;
            output << "						" << "\n" << std::flush;
            output << "									return chart;" << "\n" << std::flush;
            output << "								}" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "						});" << "\n" << std::flush;
            output << "					</script>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "			if( $Chart_Type == 'MirTail/' )" << "\n" << std::flush;
            output << "			{" << "\n" << std::flush;
            output << "				echo '<script src=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.js></script>';" << "\n" << std::flush;
            output << "				echo '<link href=http://140.113.239.177/joye/ForAgoSorting/lib/nv.d3.min.css rel=stylesheet type=text/css>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== TSV File ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$TSV = Shell_Exec( 'ls '.$Chart_Type.' | grep .tsv' );" << "\n" << std::flush;
            output << "				$TSV_List = Explode( \"\\n\", $TSV );" << "\n" << std::flush;
            output << "				$List_Size = Count( $TSV_List );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $List_Size-1; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$TSV_List[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $TSV_File == $TSV_List[$i] ) " << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					echo '>'.$TSV_List[$i].'</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceY' value='$ForceY' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== is_Abundant ====================-->" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<select name=isAbundant onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				$isAbundant_List = array('MostAbundant', 'AllmiRNA');" << "\n" << std::flush;
            output << "				$isAbundant_Size = Count( $isAbundant_List );" << "\n" << std::flush;
            output << "				if( $isAbundant == '' )" << "\n" << std::flush;
            output << "					$isAbundant = 'MostAbundant';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "				For( $i = 0; $i < $isAbundant_Size; ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					echo '<option value='.$isAbundant_List[$i].' ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					if( $isAbundant == $isAbundant_List[$i] )" << "\n" << std::flush;
            output << "						echo 'selected ';" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "					echo '>' . $isAbundant_List[$i] . '</option>';" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceY' value='$ForceY' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== miRNA Select ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$File = $Chart_Type.$TSV_File;" << "\n" << std::flush;
            output << "				$inFile = File_get_contents( $File );" << "\n" << std::flush;
            output << "				$inFile_Lines = Explode( \"\\n\", $inFile );" << "\n" << std::flush;
            output << "				$Anno_Array = Array();" << "\n" << std::flush;
            output << "				$miRNA_Array = Array();" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $j = 1; $j < Count( $inFile_Lines )-1; ++$j )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					$inFile_Line = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n" << std::flush;
            output << "					$anno = Explode( ':', $inFile_Line[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $isAbundant == 'MostAbundant' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$miRNA = Explode( '*', $anno[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( Count( $miRNA ) != 2 )" << "\n" << std::flush;
            output << "							Continue;" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					Array_Push( $Anno_Array, $anno[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $miRNA_Select == $anno[0] )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						$Tail_Array = Array();" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						for( $k = 1; $k < Count( $inFile_Line ); ++$k )" << "\n" << std::flush;
            output << "						{" << "\n" << std::flush;
            output << "							Array_Push( $Tail_Array, $inFile_Line[$k] );" << "\n" << std::flush;
            output << "						}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						$miRNA_Array[ $anno[1] ] = $Tail_Array;" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Sort( $Anno_Array, SORT_STRING );" << "\n" << std::flush;
            output << "				$Uniq_Anno_Array = Array_Unique( $Anno_Array, SORT_STRING );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<select name=miRNA_Select onchange=this.form.submit();>';" << "\n" << std::flush;
            output << "				echo '<option '; if($miRNA_Select=='') echo 'selected'; echo '>Select miRNA</option>';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < Count( $Anno_Array ); ++$i )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					if( $Uniq_Anno_Array[$i] != '' )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						echo '<option value='.$Uniq_Anno_Array[$i].' ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						if( $miRNA_Select == $Uniq_Anno_Array[$i] ) " << "\n" << std::flush;
            output << "							echo 'selected ';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "						echo '>'.$Uniq_Anno_Array[$i].'</option>';" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='ForceY' value='$ForceY' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== ForceY ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n" << std::flush;
            output << "				echo '<input type=text name=ForceY size=3 value=';" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				if( $ForceY=='' )" << "\n" << std::flush;
            output << "					echo 'Hight';" << "\n" << std::flush;
            output << "				else" << "\n" << std::flush;
            output << "					echo $ForceY;" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \" onfocus=\\\"{this.value='';}\\\">\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"</select>" << "\n" << std::flush;
            output << "					<input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='Chart_Type' value='$Chart_Type' />" << "\n" << std::flush;
            output << "					<input type='hidden' name='miRNA_Select' value='$miRNA_Select' />" << "\n" << std::flush;
            output << "					<input type='submit' value='Submit' /> " << "\n" << std::flush;
            output << "					</form>\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "#<!--================== miRNA Tail Bar Chart ====================-->" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$index = Explode( \"\\t\", $inFile_Lines[0] );" << "\n" << std::flush;
            output << "				$miRNA_Keys = Array_keys( $miRNA_Array );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				$Temp = Tempnam( '/tmp', $miRNA_Select );" << "\n" << std::flush;
            output << "				$Ftemp = fopen( $Temp, 'w' );" << "\n" << std::flush;
            output << "				Fwrite( $Ftemp, $index[0] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, \"\\t\".$miRNA_Keys[$i] );" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Fwrite( $Ftemp, \"\\n\".$index[6] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, \"\\t\".$miRNA_Array[ $miRNA_Keys[$i] ][5] );" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				Fwrite( $Ftemp, \"\\n\" );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				For( $j = 1; $j < Count( $index )-1; $j++ )" << "\n" << std::flush;
            output << "				{" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, $index[$j] );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					For( $i = 0; $i < Count( $miRNA_Keys ); $i++ )" << "\n" << std::flush;
            output << "					{" << "\n" << std::flush;
            output << "						Fwrite( $Ftemp, \"\\t\".$miRNA_Array[ $miRNA_Keys[$i] ][$j-1] );" << "\n" << std::flush;
            output << "					}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					Fwrite( $Ftemp, \"\\n\" );" << "\n" << std::flush;
            output << "				}" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				fclose( $Ftemp );" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo '<svg id=bar></svg>';" << "\n" << std::flush;
            output << "				echo \"<script>" << "\n" << std::flush;
            output << "						d3.tsv( '$Temp', function( tsv_data ) {" << "\n" << std::flush;
            output << "						var length = d3.keys( tsv_data[0] ).filter(function(key) { return key !== 'miRNA_Length'; });" << "\n" << std::flush;
            output << "						var data = new Array();" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						tsv_data.forEach( function(d) {" << "\n" << std::flush;
            output << "							d.len = length.map( function(key) {" << "\n" << std::flush;
            output << "								return {" << "\n" << std::flush;
            output << "									x: Number(key)," << "\n" << std::flush;
            output << "									y: +d[key]" << "\n" << std::flush;
            output << "								};" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "							});" << "\n" << std::flush;
            output << "							var sample = new Object();" << "\n" << std::flush;
            output << "							sample.key = d['miRNA_Length'];" << "\n" << std::flush;
            output << "							sample.values = d.len;" << "\n" << std::flush;
            output << "							data.push( sample );" << "\n" << std::flush;
            output << "						});" << "\n" << std::flush;
            output << "			" << "\n" << std::flush;
            output << "						nv.addGraph({" << "\n" << std::flush;
            output << "							generate: function() {" << "\n" << std::flush;
            output << "								var width = nv.utils.windowSize().width," << "\n" << std::flush;
            output << "									height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "								var chart = nv.models.multiBarChart()\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "					if( $ForceY != '' && $ForceY != 'Hight' )" << "\n" << std::flush;
            output << "						echo \"		.forceY([$ForceY,0])\";" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "				echo \"				.width(width)" << "\n" << std::flush;
            output << "									.height(height)" << "\n" << std::flush;
            output << "									.color( ['#000000', '#FF0000', '#0000FF', '#FFBF00', '#088A08', '#6E6E6E'])" << "\n" << std::flush;
            output << "									.stacked(true);" << "\n" << std::flush;
            output << "				" << "\n" << std::flush;
            output << "								var svg = d3.select('#bar').datum(data);" << "\n" << std::flush;
            output << "								console.log('calling chart');" << "\n" << std::flush;
            output << "								svg.transition().duration(0).call(chart);" << "\n" << std::flush;
            output << "					" << "\n" << std::flush;
            output << "								return chart;" << "\n" << std::flush;
            output << "							}," << "\n" << std::flush;
            output << "							callback: function(graph) {" << "\n" << std::flush;
            output << "									nv.utils.windowResize(function() {" << "\n" << std::flush;
            output << "										var width = nv.utils.windowSize().width;" << "\n" << std::flush;
            output << "										var height = nv.utils.windowSize().height;" << "\n" << std::flush;
            output << "										graph.width(width).height(height);" << "\n" << std::flush;
            output << "" << "\n" << std::flush;
            output << "										d3.select('#bar')" << "\n" << std::flush;
            output << "												.attr('width', width)" << "\n" << std::flush;
            output << "												.attr('height', height)" << "\n" << std::flush;
            output << "												.transition().duration(0)" << "\n" << std::flush;
            output << "												.call(graph);" << "\n" << std::flush;
            output << "									});" << "\n" << std::flush;
            output << "							}" << "\n" << std::flush;
            output << "						});" << "\n" << std::flush;
            output << "					});" << "\n" << std::flush;
            output << "				</script>\";" << "\n" << std::flush;
            output << "			}" << "\n" << std::flush;
            output << "		}" << "\n" << std::flush;
            output << "	?>" << "\n" << std::flush;
            output << "	</body>" << "\n" << std::flush;
            output << "</html>" << "\n" << std::flush;
            output.flush();
            output.close();
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
