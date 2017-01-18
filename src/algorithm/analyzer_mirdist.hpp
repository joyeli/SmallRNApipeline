#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/analyzer_ForAgoSorting.hpp>

namespace ago {
namespace algorithm {

class AnalyzerMirDist
{
    using AnaParamiRDist = AnalyzerParameter< AnalyzerTypes::miRNADistribution >;
    using AnaParamiRDistPrinter = AnalyzerParameter< AnalyzerTypes::miRDistPrinter >;
    
    using AnalyzerTypeListMir = boost::mpl::vector<
        boost::mpl::map<
              boost::mpl::pair< AnaParamiRDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::miRNADistribution >>
            , boost::mpl::pair< AnaParamiRDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParamiRDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthNameType    , boost::mpl::string< 'mi', 'RNA' >>
            , boost::mpl::pair< AnaParamiRDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParamiRDist::CalReadCountClass  , CalReadCountGMPM >
            , boost::mpl::pair< AnaParamiRDist::GetReadSeqClass    , GetReadSeed< 1, 7 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepth2NameType   , boost::mpl::string< '-1' >>
       > 
        , boost::mpl::map<
              boost::mpl::pair< AnaParamiRDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::miRNADistribution >>
            , boost::mpl::pair< AnaParamiRDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParamiRDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthNameType    , boost::mpl::string< 'mi', 'RNA' >>
            , boost::mpl::pair< AnaParamiRDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParamiRDist::GetReadSeqClass    , GetReadSeed < 1, 7 >>
            , boost::mpl::pair< AnaParamiRDist::CalReadCountClass  , CalReadCountGMOnly >
            , boost::mpl::pair< AnaParamiRDist::DbDepth2NameType   , boost::mpl::string< '-1' >>
       > 
        , boost::mpl::map<
              boost::mpl::pair< AnaParamiRDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::miRNADistribution >>
            , boost::mpl::pair< AnaParamiRDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParamiRDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParamiRDist::DbDepthNameType    , boost::mpl::string< 'mi', 'RNA' >>
            , boost::mpl::pair< AnaParamiRDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParamiRDist::GetReadSeqClass    , GetReadSeed< 1, 7 >>
            , boost::mpl::pair< AnaParamiRDist::CalReadCountClass  , CalReadCountPMOnly >
            , boost::mpl::pair< AnaParamiRDist::DbDepth2NameType   , boost::mpl::string< '-1' >>
       > 
        , boost::mpl::map<
              boost::mpl::pair< AnaParamiRDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::miRDistPrinter >>
            , boost::mpl::pair< AnaParamiRDistPrinter::AnnoIdx     , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParamiRDistPrinter::Xaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Yaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Zaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Seq, boost::mpl::string< '-2' >>>
       > 
        , boost::mpl::map<
              boost::mpl::pair< AnaParamiRDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::miRDistPrinter >>
            , boost::mpl::pair< AnaParamiRDistPrinter::AnnoIdx     , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParamiRDistPrinter::Xaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Yaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Zaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Seq, boost::mpl::string< '-2' >>>
       > 
        , boost::mpl::map<
              boost::mpl::pair< AnaParamiRDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::miRDistPrinter >>
            , boost::mpl::pair< AnaParamiRDistPrinter::AnnoIdx     , boost::mpl::int_< 2 >>
            , boost::mpl::pair< AnaParamiRDistPrinter::Xaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Yaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParamiRDistPrinter::Zaxis       , boost::mpl::pair< AnaParamiRDistPrinter::Seq, boost::mpl::string< '-2' >>>
        >
    >;

  public:

    using Analyzer = AnalyzerC< std::vector< AnnotationRawBed<> >*, AnalyzerTypeListMir >;

    Analyzer analyzer;

    AnalyzerMirDist()
    {}

    void tailing_ratio( std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result )
    {
        std::map< std::string, std::map< std::string, double >> TR_table;
        std::map< std::string, std::map< std::string, double >>::iterator tableIt;
        std::map< std::string, double >::iterator lenIt;

        for( auto& result_table : Analyzer_Result )
        {
            tableIt = result_table.find( ".MirDist_GMPM_read_count" );

            if( tableIt != result_table.end() )
            {
                int flag = 0;
                for( auto& table : result_table )
                {
                    if( flag == 0 )
                        TR_table.emplace( ".MirDist_Tailing_Ratio", table.second );
                    else
                        TR_table.emplace( table.first, table.second );

                    flag++;
                }
            }
        }

        for( auto& result_table : Analyzer_Result )
        {
            tableIt = result_table.find( ".MirDist_PM_read_count" );

            if( tableIt != result_table.end() )
            {
                int flag = 0;
                for( auto& TR : TR_table )
                {
                    if( flag != 0 )
                    {
                        tableIt = result_table.find( TR.first );

                        if( tableIt != result_table.end() )
                        {
                            for( auto& len : TR.second )
                            {
                                lenIt = tableIt->second.find( len.first );

                                if( lenIt != tableIt->second.end() )
                                {
                                    if( len.second != 0 )
                                        len.second = lenIt->second / len.second * 100;
                                }
                            }
                        }
                    }

                    flag++;
                }
            }
        }

        Analyzer_Result.push_back( TR_table );
    }

    void tailing_ratio_for_PM( std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result,
            std::vector< AnnotationRawBed<> >& Annotation_RawBed,
            std::map< std::string, std::string >& Genome_Table )
    {
        std::map< std::string, std::map< std::string, double >> TR_table;

        std::map< int, std::map< char, double >> tailing_map;
        std::map< int, std::map< char, double >>::iterator tailing_It;

        for( auto& AnnoRawBed : Annotation_RawBed )
        {
            if( AnnoRawBed.is_filtered_ == 0 && AnnoRawBed.getTail().size() != 0 )
            {
                for( auto& AnnoRawBed_Info : AnnoRawBed.annotation_info_ )
                {
                    for( int i = 0; i < AnnoRawBed_Info.size(); i+=2 )
                    {
                        if( AnnoRawBed_Info[i] == "miRNA" )
                        {
                            std::map< char, double > tail_map;

                            int length = AnnoRawBed.length_ - AnnoRawBed.tail_length_;
                            double read_count = AnnoRawBed.reads_count_ / AnnoRawBed.multiple_alignment_site_count_;
                            std::string tail = AnnoRawBed.getTail();

                            tailing_It = tailing_map.find( length );

                            char N( tail[0] );

                            for( int i = 0; i < tail.length(); ++i )
                            {
                                if( N != tail[i] )
                                {
                                    N = 'X';
                                    break;
                                }
                            }

                            if( tailing_It != tailing_map.end() )
                            {
                                if( tailing_It->second.find( N ) != tailing_It->second.end() )
                                    tailing_It->second.find( N )->second = read_count + tailing_It->second.find( N )->second;
                                else
                                    tailing_It->second.emplace( N, read_count );
                            }
                            else
                            {
                                tail_map.emplace( 'A', 0 );
                                tail_map.emplace( 'C', 0 );
                                tail_map.emplace( 'G', 0 );
                                tail_map.emplace( 'T', 0 );
                                tail_map.emplace( 'X', 0 );
                                tail_map.emplace( N, read_count );
                                tail_map.find(N)->second = read_count + tail_map.find(N)->second;
                                tailing_map.emplace( length, tail_map );
                            }
                        }
                    }
                }
            }
        }

        for( auto& tail_length : tailing_map )
        {
            double total(0);

            for( auto& tail_N : tail_length.second )
                total = tail_N.second + total;

            for( auto& tail_N : tail_length.second )
                tail_N.second = tail_N.second / total;
        }

        std::string NN;
        std::map< std::string, double > table;
        TR_table.emplace( ".PM_Tailing_Ratio", table );

        for( auto& tail_length : tailing_map )
        {
            for( auto& tail_N : tail_length.second )
            {
                NN.push_back( tail_N.first );
                table.emplace( NN, tail_N.second );
                NN.clear();
            }

            TR_table.emplace( std::to_string( tail_length.first ), table );
            table.clear();
        }

        Analyzer_Result.push_back( TR_table );
    }

    void tailing_ratio_for_miRNA( std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result,
            std::vector< AnnotationRawBed<> >& Annotation_RawBed,
            std::map< std::string, std::string >& Genome_Table )
    {
        std::map< std::string, std::map< std::string, double >> TR_table;

        std::map< std::string, std::map< int, std::map< char, double >>> tailing_map;
        std::map< std::string, std::map< int, std::map< char, double >>>::iterator anno_It;

        for( auto& AnnoRawBed : Annotation_RawBed )
        {
            if( AnnoRawBed.is_filtered_ == 0 )
            {
                for( auto& AnnoRawBed_Info : AnnoRawBed.annotation_info_ )
                {
                    for( int i = 0; i < AnnoRawBed_Info.size(); i+=2 )
                    {
                        if( AnnoRawBed_Info[i] == "miRNA" )
                        {
                            std::string anno( AnnoRawBed_Info[i+1] );
                            std::string sequ( AnnoRawBed.getReadSeq( Genome_Table ));
                            std::string anno_seed( anno + "_" + sequ.substr( 1, 7 ));

                            std::string tail = AnnoRawBed.getTail();
                            int length = AnnoRawBed.length_ - AnnoRawBed.tail_length_;

                            double read_count = AnnoRawBed.reads_count_ / AnnoRawBed.multiple_alignment_site_count_;

                            std::map< int, std::map< char, double >> len_map;
                            std::map< int, std::map< char, double >>::iterator len_It;
                            std::map< char, double > tail_map;

                            char N( tail[0] );

                            for( int i = 0; i < tail.length(); ++i )
                            {
                                if( N != tail[i] )
                                {
                                    N = 'X';
                                    break;
                                }
                            }

                            if( tail.size() == 0 )
                                N = 'Z';

                            anno_It = tailing_map.find( anno_seed );

                            if( anno_It != tailing_map.end() )
                            {
                                len_It = anno_It->second.find( length );

                                if( len_It != anno_It->second.end() )
                                {
                                    if( len_It->second.find( N ) != len_It->second.end() )
                                    {
                                        len_It->second.find( N )->second = read_count + len_It->second.find( N )->second;
                                    }
                                    else
                                    {
                                        len_It->second.emplace( N, read_count );
                                    }
                                }
                                else
                                {
                                    tail_map.emplace( 'A', 0 );
                                    tail_map.emplace( 'C', 0 );
                                    tail_map.emplace( 'G', 0 );
                                    tail_map.emplace( 'T', 0 );
                                    tail_map.emplace( 'X', 0 );
                                    tail_map.emplace( 'Z', 0 );
                                    tail_map.find(N)->second = read_count + tail_map.find(N)->second;
                                    anno_It->second.emplace( length, tail_map );
                                }
                            }
                            else
                            {
                                tail_map.emplace( 'A', 0 );
                                tail_map.emplace( 'C', 0 );
                                tail_map.emplace( 'G', 0 );
                                tail_map.emplace( 'T', 0 );
                                tail_map.emplace( 'X', 0 );
                                tail_map.emplace( 'Z', 0 );
                                tail_map.find(N)->second = read_count + tail_map.find(N)->second;
                                len_map.emplace( length, tail_map );
                                tailing_map.emplace( anno_seed, len_map );
                            }
                        }
                    }
                }
            }
        }

        for( auto& anno_seed : tailing_map )
        {
            for( auto& length : anno_seed.second )
            {
                double total(0);

                for( auto& tail_N : length.second )
                    total = tail_N.second + total;

                if( total == 0 )
                    continue;

                for( auto& tail_N : length.second )
                    tail_N.second = tail_N.second / total;
            }
        }

        std::string NN;
        std::map< std::string, double > table;
        TR_table.emplace( ".miRNA_Tailing_Ratio", table );

        for( auto& anno_seed : tailing_map )
        {
            for( auto& length : anno_seed.second )
            {
                for( auto& tail_N : length.second )
                {
                    NN.push_back( tail_N.first );
                    table.emplace( NN, tail_N.second );
                    NN.clear();
                }

                TR_table.emplace( anno_seed.first + ":" + std::to_string( length.first ), table );
                table.clear();
            }
        }

        Analyzer_Result.push_back( TR_table );
    }

    void tailing_detail_for_miRNA(
              std::pair< std::string, std::vector< AnnotationRawBed<> >>& sample_rawbed
            , std::map< std::string, std::string >& genome_table
            , std::string& output_path
    )
    {
        std::map< std::string, std::map< size_t, std::map< std::string, double >>> anno_tail_map;
        //        anno_seed               lens              tail_seq    count

        std::map< std::string, std::map< size_t, std::map< std::string, double >>>::iterator anno_tail_it;
        std::map< size_t, std::map< std::string, double >>::iterator len_tail_it;
        std::map< std::string, double >::iterator tail_it;

        for( auto& raw_bed : sample_rawbed.second )
        {
            if( raw_bed.is_filtered_ == 0 )
            {
                for( auto& raw_bed_info : raw_bed.annotation_info_ )
                {
                    for( int i = 0; i < raw_bed_info.size(); i+=2 )
                    {
                        if( raw_bed_info[i] == "miRNA" )
                        {
                            std::string anno( raw_bed_info[i+1] );
                            std::string sequ( raw_bed.getReadSeq( genome_table ));
                            std::string anno_seed( anno + "_" + sequ.substr( 1, 7 ));

                            std::string tail = raw_bed.getTail();
                            size_t length = raw_bed.length_ - raw_bed.tail_length_;

                            double read_count = raw_bed.reads_count_ / raw_bed.multiple_alignment_site_count_;
                            
                            anno_tail_it = anno_tail_map.find( anno_seed );

                            if( anno_tail_it != anno_tail_map.end() )
                            {
                                len_tail_it = anno_tail_it->second.find( length );

                                if( len_tail_it != anno_tail_it->second.end() )
                                {
                                    tail_it = len_tail_it->second.find( tail );

                                    if( tail_it != len_tail_it->second.end() )
                                    {
                                        tail_it->second += read_count;
                                    }
                                    else
                                    {
                                        len_tail_it->second.emplace( tail, read_count );
                                    }
                                }
                                else
                                {
                                    std::map< std::string, double > tail_temp;
                                    tail_temp.emplace( tail, read_count );

                                    anno_tail_it->second.emplace( length, tail_temp );
                                }
                            }
                            else
                            {
                                std::map< std::string, double > tail_temp;
                                tail_temp.emplace( tail, read_count );

                                std::map< size_t, std::map< std::string, double >> len_tail_temp;
                                len_tail_temp.emplace( length, tail_temp );

                                anno_tail_map.emplace( anno_seed, len_tail_temp );
                            }
                        }
                    }
                }
            }
        }

        std::ofstream output( output_path + "/" + sample_rawbed.first + "_tailing.tsv" );
        output << "AnnoSeed\tLen\tTailSeq\tCount\n";

        for( auto& anno_len : anno_tail_map )
        {
            for( auto& len_tail : anno_len.second )
            {
                for( auto& tail_count : len_tail.second )
                {
                    output
                        << anno_len.first << "\t"
                        << len_tail.first << "\t"
                        << tail_count.first << "\t"
                        << tail_count.second << "\n"
                        ;
                }
            }
        }

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
