#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <AGO/algorithm/analyzer_ForAgoSorting.hpp>

namespace ago {
namespace algorithm {

class AnalyzerLenDist
{
    using AnaParaLenDist = AnalyzerParameter< AnalyzerTypes::LengthDistribution >;
    using AnaParaLenDistPrinter = AnalyzerParameter< AnalyzerTypes::LenDistPrinter >;
    
    using AnalyzerTypeListLen = boost::mpl::vector<
        boost::mpl::map<
              boost::mpl::pair< AnaParaLenDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::LengthDistribution >>
            , boost::mpl::pair< AnaParaLenDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParaLenDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthNameType    , boost::mpl::string< '-1' >>
            , boost::mpl::pair< AnaParaLenDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParaLenDist::CalReadCountClass  , CalReadCountGMPM >
        >
        , boost::mpl::map<
              boost::mpl::pair< AnaParaLenDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::LengthDistribution >>
            , boost::mpl::pair< AnaParaLenDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParaLenDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthNameType    , boost::mpl::string< '-1' >>
            , boost::mpl::pair< AnaParaLenDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParaLenDist::CalReadCountClass  , CalReadCountGMOnly >
        >
        , boost::mpl::map<
              boost::mpl::pair< AnaParaLenDist::AnalyzerType       , boost::mpl::int_< AnalyzerTypes::LengthDistribution >>
            , boost::mpl::pair< AnaParaLenDist::FilterType         , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParaLenDist::DbIndexType        , boost::mpl::int_< -1 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthType        , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParaLenDist::DbDepthNameType    , boost::mpl::string< '-1' >>
            , boost::mpl::pair< AnaParaLenDist::GetReadLengthClass , GetReadPrefixLength >
            , boost::mpl::pair< AnaParaLenDist::CalReadCountClass  , CalReadCountPMOnly >
        >
        , boost::mpl::map<
              boost::mpl::pair< AnaParaLenDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::LenDistPrinter >>
            , boost::mpl::pair< AnaParaLenDistPrinter::AnnoIdx     , boost::mpl::int_< 0 >>
            , boost::mpl::pair< AnaParaLenDistPrinter::Xaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Yaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Zaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Seq, boost::mpl::string< '-2' >>>
        >
        , boost::mpl::map<
              boost::mpl::pair< AnaParaLenDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::LenDistPrinter >>
            , boost::mpl::pair< AnaParaLenDistPrinter::AnnoIdx     , boost::mpl::int_< 1 >>
            , boost::mpl::pair< AnaParaLenDistPrinter::Xaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Yaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Zaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Seq, boost::mpl::string< '-2' >>>
        >
        , boost::mpl::map<
              boost::mpl::pair< AnaParaLenDistPrinter::AnalyzerType, boost::mpl::int_< AnalyzerTypes::LenDistPrinter >>
            , boost::mpl::pair< AnaParaLenDistPrinter::AnnoIdx     , boost::mpl::int_< 2 >>
            , boost::mpl::pair< AnaParaLenDistPrinter::Xaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Anno, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Yaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Len, boost::mpl::string< '-1' >>>
            , boost::mpl::pair< AnaParaLenDistPrinter::Zaxis       , boost::mpl::pair< AnaParaLenDistPrinter::Seq, boost::mpl::string< '-2' >>>
        >
    >;

  public:

	using Analyzer = AnalyzerC< std::vector< AnnotationRawBed<> >*, AnalyzerTypeListLen >;

    Analyzer analyzer;

    AnalyzerLenDist()
    {}

	void tailing_ratio( std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result )
	{
		std::map< std::string, std::map< std::string, double >> TR_table;
		std::map< std::string, std::map< std::string, double >>::iterator tableIt;
		std::map< std::string, double >::iterator lenIt;

		for( auto& result_table : Analyzer_Result )
		{
			tableIt = result_table.find( ".LenDist_GMPM_read_count" );

			if( tableIt != result_table.end() )
			{
				int flag = 0;
				for( auto& table : result_table )
				{
					if( flag == 0 )
						TR_table.emplace( ".LenDist_Tailing_Ratio", table.second );
					else
						TR_table.emplace( table.first, table.second );

					flag++;
				}
			}
		}

		for( auto& result_table : Analyzer_Result )
		{
			tableIt = result_table.find( ".LenDist_PM_read_count" );

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
};

} // end of namespace algorithm
} // end of namespace ago
