#ifndef MIRDISPRINTER_HPP_
#define MIRDISPRINTER_HPP_

#include <fstream>
#include <array>

#include <boost/filesystem.hpp>


#include <pokemon/constant_def.hpp>
#include <pokemon/tuple_utility.hpp>
#include <pokemon/iohandler/iohandler.hpp>
#include <pokemon/iohandler/basespace_def.hpp>
#include "boost/algorithm/string/split.hpp"

template<>
class AnalyzerParameter<AnalyzerTypes::miRDistPrinter>
{
public:
	/** 
	 * @brief 定義 length distribution參數有哪些
	 */
	
	/// @brief AnalyzerType 特化Analyzer的type，每一個AnalyzerParameter都應該要有此參數
	typedef boost::mpl::int_<0> AnalyzerType;
	
	/// @brief 
	typedef boost::mpl::int_<1> AnnoIdx;
	
	typedef boost::mpl::int_<2> Xaxis;
	
	typedef boost::mpl::int_<3> Yaxis;
	
	typedef boost::mpl::int_<4> Zaxis;
	
	typedef boost::mpl::int_<5> Xlimit;
	
	typedef boost::mpl::int_<6> Ylimit;
	
	typedef boost::mpl::int_<7> Zlimit;
	
	typedef boost::mpl::int_<8> Len;
	
	typedef boost::mpl::int_<9> Anno;
	
	typedef boost::mpl::int_<10> Seq;
	
	typedef boost::mpl::int_<11> PrefixName;
};


template<class INPUT_TYPE, class ANALYZER_TYPELIST>
class AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRDistPrinter>
{
public:
	/// @brief 此為要使用的 analyzer parameter
	typedef AnalyzerParameter <AnalyzerTypes::miRDistPrinter> AnaPara;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::AnnoIdx, boost::mpl::int_<0> >::type AnnoIdx;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Xaxis, boost::mpl::pair< AnaPara::Len, boost::mpl::string<'-1'> > >::type Xaxis;
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Yaxis, boost::mpl::pair< AnaPara::Anno, boost::mpl::string<'-1'> > >::type Yaxis;
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Zaxis, boost::mpl::pair< AnaPara::Seq, boost::mpl::string<'-2'> > >::type Zaxis;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Xlimit, boost::mpl::int_<0> >::type Xlimit;
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Ylimit, boost::mpl::int_<0> >::type Ylimit;
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Zlimit, boost::mpl::int_<0> >::type Zlimit;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::PrefixName, boost::mpl::string<'Tab', 'le'> >::type PrefixName;
	
	typedef std::vector< std::map<int, std::map<std::tuple<int,int,int,int,std::string,std::string>, std::tuple<double, double, double> > > > gOutSetType ;
	gOutSetType &gOutSet;
	
	std::vector<std::string> &gOutNamePrefix;
	
	std::string output_path;
	
	/// @brief 輸入的資料
	INPUT_TYPE &in;
	std::string GMPM; //Joye
	
	// Z(X(Y))
	std::map< std::string, std::map< std::pair< std::string, std::string >, std::map< std::string, double >>> cal_table;//Joye
//Joye	std::map< std::string, std::map< std::string, std::map< std::string, double >>> cal_table;

	/// @brief tuple< filter = -1, filter = 1, read_count >
	std::map< std::string, std::vector< double >> cal_normalize;
	
	/// @memberof AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>
	/// @brief AnalyzerImpl 建構子
	AnalyzerImpl(INPUT_TYPE &i)
		: in( i )
		, gOutSet( AnalyzerImpl< INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution >::gOutSet_ )
		, gOutNamePrefix( AnalyzerImpl< INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution >::gOutNamePrefix_ )
		, output_path("")
	{}
	
	void* operator()( int this_analyzer_count, size_t pipe_index, bool eof_flag, std::string& OutputPath, std::string& SampleName, std::map< std::string, std::string >& Genome_Table, std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result, int barcode_index=0 )
	{
		if( eof_flag )
		{
			std::map< std::string, std::map< std::string, double >> LenMir_Map; //Joye

			output_path = OutputPath + "/" + SampleName + "/"; //Joye
			
			Analysis();
			mkdir();
			
			{
				std::map< std::string, std::map< std::string, double >> LenMir_Map; //Joye

				output( LenMir_Map, "read_count" ); //Joye
				Analyzer_Result.push_back( LenMir_Map ); //Joye

			}
			{
				std::map< std::string, std::map< std::string, double >> LenMir_Map; //Joye

				output( LenMir_Map, "ppm", 1 ); //Joye
				Analyzer_Result.push_back( LenMir_Map ); //Joye
			}

//Joye			output( "read_count" );
//Joye			output( "ppm", 1 );
			
			//std::cout << "filter seqdepth " << std::get<0>(gOutSet[0][-1][std::make_tuple(0,1,0,0,"","")]) << std::endl;
			//std::cout << "All seqdepth " << std::get<0>(gOutSet[0][-1][std::make_tuple(0,-1,0,0,"","")]) << std::endl;			
		}

		if( AnnoIdx::value+1 == gOutSet.size() )	//Joye
		{
			AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutSet_.clear();	//Joye
			AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutNamePrefix_.clear();	//Joye
		}

		return (void*) &in;
	}

	std::string get_XYZ_name( const int i )
	{
		if( i == 8)
			return "LEN";
		else if( i == 9 )
			return "ANNO";
		else if( i == 10 )
			return "SEQ";
	}

	std::string get_filename( const std::string& Z_name, const std::string tail_name )
	{
		std::string NamePrefix = gOutNamePrefix[AnnoIdx::value];

		std::vector< std::string > tmp_vv; //Joye
		boost::split( tmp_vv, NamePrefix, boost::is_any_of( "." )); //Joye
		GMPM = tmp_vv[3]; //Joye

/*Joye
		if (Z_name[0] == '0')
		{
			std::vector<std::string> tmp_v;
			boost::split( tmp_v, NamePrefix, boost::is_any_of( "." ));
			NamePrefix = "";
			for(int i=0; i != tmp_v.size()-2; i++)
			{
				NamePrefix += tmp_v[i]+".";
			}
			NamePrefix += "merged.no_detail";
		}
*/

		std::string filename(boost::mpl::c_str<PrefixName>::value);
		filename += ".";
		filename += NamePrefix + "."; //lendist.Normal.full.GMPM.biotype.no_detail
//Joye		filename += Z_name + "_";//2.1.0.0
//Joye		filename += get_XYZ_name(Zaxis::first::value) + ".";//Total_SEQ.
		filename += get_XYZ_name(Xaxis::first::value) + "vs";//LENvsANNO.
		filename += get_XYZ_name(Yaxis::first::value) + ".";//read_count.
		filename += tail_name;
		return filename; 
	}

	void Analysis()
	{
		auto &OutSet = gOutSet[AnnoIdx::value];

		for( auto &len_pair : OutSet )
		{
			int len = len_pair.first;

			if(len == -1)
				continue;

			for( auto &anno_value_pair : len_pair.second )
			{
				const std::string &anno = std::get<4>( anno_value_pair.first );
				const std::string &seq = std::get<5>( anno_value_pair.first );
				
				double &read_count = std::get<0>(anno_value_pair.second);
				
				int dim_x = Xaxis::first::value-8;
				int dim_y = Yaxis::first::value-8;
				int dim_z = Zaxis::first::value-8;
				
				std::string type_x( boost::mpl::c_str< typename Xaxis::second >::value );
				std::string type_y( boost::mpl::c_str< typename Yaxis::second >::value );
				std::string type_z( boost::mpl::c_str< typename Zaxis::second >::value );
				
				std::array< std::string, 3 > origin_type({{ type_x, type_y, type_z }});
				
				std::string X = get_key( origin_type[dim_x], std::to_string(len) );
				std::string Y = get_key( origin_type[dim_y], anno );
				std::string Z = get_key( origin_type[dim_z], seq );
				std::string S = get_key( origin_type[dim_x], seq );

				std::string sum_x_name = "SUM_" + get_XYZ_name( Xaxis::first::value );
				std::string sum_y_name = "SUM_" + get_XYZ_name( Yaxis::first::value );
				
				std::array< std::string, 4 > XYZ ({{ X, Y, Z, S }});//Joye
//Joye				std::array< std::string, 3 > XYZ ({{ X, Y, Z }});
				
				std::string prefix("");
				prefix += std::to_string( std::get<0>( anno_value_pair.first ) ) + "."; // sys
				prefix += std::to_string( std::get<1>( anno_value_pair.first ) ) + "."; // filter
				prefix += std::to_string( std::get<2>( anno_value_pair.first ) ) + "."; // db_idx
				prefix += std::to_string( std::get<3>( anno_value_pair.first ) ) + "."; // db_depth

				cal_table[ prefix + XYZ[dim_z] ][ std::make_pair( XYZ[dim_x], XYZ[3] )][ XYZ[dim_y] ] += read_count;//Joye
				cal_table[ prefix + XYZ[dim_z] ][ std::make_pair( sum_x_name, sum_x_name )][ sum_y_name ] += read_count;//Joye
				cal_table[ prefix + XYZ[dim_z] ][ std::make_pair( XYZ[dim_x], XYZ[3] )][ sum_y_name ] += read_count;//Joye
				cal_table[ prefix + XYZ[dim_z] ][ std::make_pair( sum_x_name, sum_x_name )  ][ XYZ[dim_y] ] += read_count;//Joye

/*Joye
				cal_table[ prefix + XYZ[dim_z] ][ XYZ[dim_x] ][ XYZ[dim_y] ] += read_count;
				cal_table[ prefix + XYZ[dim_z] ][ sum_x_name ][ sum_y_name ] += read_count;
				cal_table[ prefix + XYZ[dim_z] ][ XYZ[dim_x] ][ sum_y_name ] += read_count;
				cal_table[ prefix + XYZ[dim_z] ][ sum_x_name  ][ XYZ[dim_y] ] += read_count;
*/	
		
				//record ppm
				if( cal_normalize.find( prefix + XYZ[dim_z] ) == cal_normalize.end() )
				{
					int db_idx = std::get<2>( anno_value_pair.first );
					double normalize_filter = std::get<0>( gOutSet[0][-1][ std::make_tuple( 0,1,db_idx,0,"","" ) ]); // 1
					double normalize_all = std::get<0>( gOutSet[0][-1][ std::make_tuple(0,-1,db_idx,0,"","" ) ]); // -1
					cal_normalize.insert({ prefix + XYZ[dim_z], { normalize_all, normalize_filter, 1000000 } });
				}
			}
		}
	}

	void mkdir()
	{
		boost::filesystem::path dir_output("output");
		boost::filesystem::path dir_analyzer(output_path);
		boost::filesystem::create_directory(dir_analyzer);
	}
	
	void output( std::map< std::string, std::map< std::string, double >>& LenMir_Map, const std::string &unit = "read_count", double normalize_idx = 2 )
	{
		std::vector< std::shared_ptr< IoHandlerBaseSpace >> vec;

		int count_z(0);
		for(auto &dim_z : cal_table)
		{
			double normalize = cal_normalize[dim_z.first][normalize_idx];
			
			if(count_z > Zlimit::value && Zlimit::value != 0)
				break;

			std::ofstream* out_ = new std::ofstream (output_path + get_filename(dim_z.first, unit + ".tsv"));
//Joye			std::cout << output_path + get_filename(dim_z.first, unit + ".tsv") << std::endl;
			*(out_) << "====================" << std::endl;
			*(out_) << get_filename(dim_z.first, unit + ".tsv") << std::endl;
			*(out_) << "--------------------" << std::endl;
			
			std::vector< std::pair< std::string, double> > all_title;
			int count_x(0);

			for(auto &dim_x : dim_z.second)
			{
				if(count_x > Xlimit::value && Xlimit::value != 0)
					break;
				
				if(count_x == 0)
				{
					std::string LenMir_Name( ".MirDist_" + GMPM + "_" + unit ); //Joye
					std::map< std::string, double > LenMir_Length; //Joye

					*(out_) << "ANNO" << "\t";//Joye
					*(out_) << "SEED" << "\t";//Joye
//Joye					*(out_) << "table" << "\t";

					std::string sum_x_name = "SUM_" + get_XYZ_name(Xaxis::first::value);	

					for(auto &dim_y : (dim_z.second)[std::make_pair( sum_x_name, sum_x_name )])//Joye
//Joye					for(auto &dim_y : (dim_z.second)[sum_x_name])
					{
						all_title.push_back(dim_y);

						LenMir_Length.emplace( dim_y.first, atof(dim_y.first.c_str()) ); //Joye
					}// 將y軸推入all_title

					LenMir_Map.emplace( LenMir_Name, LenMir_Length ); //Joye

					// sort y (anno) big > small
					std::sort(all_title.begin(), all_title.end(), 
						[](const std::pair<std::string, double> &a, const std::pair<std::string, double> &b)
						{
							return a.second > b.second;
						}
					);//排序y軸有小到大

					for(auto &title : all_title)
					{
						if(title.first == "")
							*(out_) << "NO_SEQ" << "\t";
						else
							*(out_) << title.first << "\t";//印出title
					}
						
					*(out_) << "\n" ;
				}

				std::string LenMir_Name( dim_x.first.first + "_" + dim_x.first.second ); //Joye
				std::map< std::string, double > LenMir_Length; //Joye

				*(out_) << dim_x.first.first << "\t";//Joye
				*(out_) << dim_x.first.second << "\t";//Joye
//Joye				*(out_) << dim_x.first << "\t";
				
				int count_y(0);
				for(auto &title : all_title)
				{
					if(count_y > Ylimit::value && Ylimit::value != 0)
						break;
					out_->precision(3);
					*(out_) << std::fixed << ( (dim_x.second)[title.first] * 1000000 / normalize ) << "\t";
					++count_y;

					LenMir_Length.emplace( title.first, ( (dim_x.second)[title.first] * 1000000 / normalize ) ); //Joye
				}

				LenMir_Map.emplace( LenMir_Name, LenMir_Length ); //Joye

				*(out_) << "\n" ;
				++count_x;
			}
			*(out_) << "--------------------" << std::endl;
			++count_z;
			out_->close ();
		}
	}
	
	std::string get_key( const std::string &type, const std::string &value )
	{
		if(type == "-1")
		{
			return value;
		}
		else if(type == "-2")
		{
			return "Total";
		}
		else if(type == "-3")
		{
			if(value == "")
				return "No-seq";
			else
				return "Having-seq";
		}
		else
		{
			if(type == value)
				return value;
			else
				return "Other";
		}
	}
};
#endif
