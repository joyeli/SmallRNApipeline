#ifndef ANALYZER_MIRNADISTRIBUTION_HPP_
#define ANALYZER_MIRNADISTRIBUTION_HPP_

#include <cassert>

#include <pokemon/constant_def.hpp>
#include <pokemon/tuple_utility.hpp>

#include <pokemon/analyzer/analyzer_utility.hpp>
#include <AGO/algorithm/analyzer_policy_ForAgoSorting.hpp>
//#include "analyzer_printer_old.hpp"

/**
 * @brief Analyzer 此class為提供計算Langth distribution的Analyzer參數
 * @tparam ANALYZER_TYPE 特化為 miRNADistribution
 */
template<>
class AnalyzerParameter<AnalyzerTypes::miRNADistribution>
{
public:
	/** 
	 * @brief 定義 length distribution參數有哪些
	 */
	
	/// @brief AnalyzerType 特化Analyzer的type，每一個AnalyzerParameter都應該要有此參數
	typedef boost::mpl::int_<0> AnalyzerType;
	
	/// @brief FilterType 決定Filter後，此analyzer 要取那個Tag，-1=>全，1=>去掉filter tag=1，0=>卻掉filter tag=0。每一個AnalyzerParameter都應該要有此參數
	typedef boost::mpl::int_<1> FilterType;
	
	/// @brief DbIndexType 決定annotation 的 db。可以為空，代表不做。可以為 -1，代表全做。可以為數字，指定db。
	typedef boost::mpl::int_<2> DbIndexType;
	
	/// @brief DbDepthType 決定 第N個 annotation。可以為空，代表不做。可以為 -1，代表全做。可以為數字，指定第N個 annotation。
	typedef boost::mpl::int_<3> DbDepthType;
	
	/// @brief DbDepthNameType 決定annotation name。可以為空，代表全部名字。可以為字串(boost::mpl::string)，指定 annotation name為何。
	typedef boost::mpl::int_<4> DbDepthNameType;
	
	/// @brief GetReadLengthClass 決定取得 anno raw bed length 的方法，也就是length distribution 的 length，預設回傳 read length，也可以自定回傳其他數值
	typedef boost::mpl::int_<5> GetReadLengthClass;
	
	/// @brief CalReadCountClass 決定取得 anno raw bed count 的方法，也就是length distribution 的 value，預設回傳 read count，也可以自定回傳其他數值
	typedef boost::mpl::int_<6> CalReadCountClass;
	
	/// @brief GetReadSeqClass 決定取得 anno raw bed sequence 的方法，也就是read sequence 的 content，預設回傳 空字串，也可以自定回傳其他數值
	typedef boost::mpl::int_<7> GetReadSeqClass;
	
	typedef boost::mpl::int_<8> Printer;

	typedef boost::mpl::int_<-1> ReturnType;
	
	typedef boost::mpl::int_<9> DbDepth2NameType;
};


class ANALYZER_TYPELIST_GLOBAL_MIRDIST {};


template<class INPUT_TYPE>
class AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>
{
public:
	typedef std::map<int, std::map<std::tuple<int,int,int,int,std::string,std::string>, std::tuple<double, double, double> > > OutPutType ;
	
	/// @brief 因為每個run為mt，所以最後要整合進一個 global 的資料結構
	static std::vector<OutPutType> gOutSet_;
	static std::vector<std::string> gOutNamePrefix_;
	
	static std::mutex gOutMutex_;
	
	static void ClearContent (void)
	{
		gOutSet_.clear();
	}
};


/**
 * @brief Analyzer 的實作，此為計算 Length distribation。幾乎範型的可以計算所有 length distribution，可以用參數設定要的長度分布
 * @tparam INPUT_TYPE 輸入資料的型別，一定為 vector，通常為 vector<AnnoRawBed>
 * @tparam ANALYZER_TYPELIST Analyzer要用的參數設定，通常為 boost::mpl::map<boost::mpl::pair<KEY, VALUE> >
 * @tparam ANALYZER_TYPE 特化 Analyzer用的參數，此為 Length distribation 特化
 */
template<class INPUT_TYPE, class ANALYZER_TYPELIST>
class AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>
{
public:
	/// @brief 此為要使用的 analyzer parameter
	typedef AnalyzerParameter <AnalyzerTypes::miRNADistribution> AnaPara;
	
	/// @brief FilterType 決定Filter後，此analyzer 要取那個Tag，-1=>全，1=>去掉"filter tag"=1，0=>去掉 "filter tag"=0。每一個AnalyzerParameter都應該要有此參數
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::FilterType, boost::mpl::int_<1> >::type FilterType;
	
	/// @brief DbIndexType 決定annotation 的 db。可以為空，代表不做，此轉為 -2。可以為 -1，代表全做。可以為數字，指定db。
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::DbIndexType, boost::mpl::int_<-2> >::type DbIndexType;
	
	/// @brief DbDepthType 決定 第N個 annotation。可以為空，代表不做，此轉為 -2。可以為 -1，代表全做。可以為數字，指定第N個 annotation。
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::DbDepthType, boost::mpl::int_<-2> >::type DbDepthType;
	
	/// @brief DbDepthNameType 決定annotation name。可以為空，此轉為 "-1"，代表全部名字。可以為字串(boost::mpl::string)，指定 annotation name為何。
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::DbDepthNameType, boost::mpl::string<'-1'> >::type DbDepthNameType;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::DbDepth2NameType, boost::mpl::string<'-3'> >::type DbDepth2NameType;
	
	/// @brief GetReadLengthClass 決定取得 anno raw bed length 的方法，也就是length distribution 的length，預設回傳 read length，也可以自定回傳其他數值
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::GetReadLengthClass, GetReadLengthDefault >::type GetReadLengthClass;
	
	/// @brief CalReadCountClass 決定取得 anno raw bed count 的方法，也就是length distribution 的 value，預設回傳 read count，也可以自定回傳其他數值
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::CalReadCountClass, CalReadCountDefault >::type CalReadCountClass;
	
	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::GetReadSeqClass, GetReadSeqDefault >::type GetReadSeqClass;
	
	//typedef typename at<ANALYZER_TYPELIST, typename AnaPara::Printer, LenDistSeqPrinter >::type PrinterClass;

	typedef typename at<ANALYZER_TYPELIST, typename AnaPara::ReturnType >::type ReturnType;
	
	/// @brief 記錄此run為第N個使用者要求的length distribution，因為輸出資料節過為 vector<map<...>>
	int this_analyzer_count_;
	
	/// @brief 輸入的資料
	INPUT_TYPE &in;
	
	/** 
	 * @brief 輸出的資料型別
	 * \n output vector {read_length, tuple{sys, filter, db_idx, db_depth, db_depth_name, ReadSeqinfo}, tuple< read_number, all ppm, filter ppm >} 
	 * \n sys 此為特殊參數，主要是output所必須記錄，0=> db_index為空 and db_depth為空，1=> db_index為空 and db_depth不為空，2=>db_index不為空 and db_depth不為空
	 * \n sys 是因為db_index, db_depth "空"沒辦法記錄。ps. 自從空改成 -2後，應該就可以不需要 sys了，但目前還有待修改...
	 */
	typedef std::map< int, std::map< std::tuple< int,int,int,int,std::string,std::string >, std::tuple< double, double, double >>> OutPutType ;
	/// @brief 輸出的資料
	OutPutType out_set_;
	
	/// @brief 因為每個run為mt，所以最後要整合進一個 global 的資料結構
	static std::vector< OutPutType > &gOutSet_;
	static std::vector<std::string> &gOutNamePrefix_;
	
	static std::mutex &gOutMutex_;
	
	/// @memberof AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>
	/// @brief AnalyzerImpl 建構子
	AnalyzerImpl(INPUT_TYPE &i)
		: in (i), this_analyzer_count_(0)
	{}
	/// @memberof AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>
	/// @brief AnalyzerImpl 建構子
	AnalyzerImpl()
	{}
	
	/**
	 * @memberof AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>
	 * @brief 主要執行Analyzer的對外介面
	 * @param this_analyzer_count 第N個 使用者要求的length distribution，因為輸出資料節過為 vector<map<...>>
	 * @return void
	 */
	void* operator()( int this_analyzer_count, size_t pipe_index, bool eof_flag, std::string& OutputPath, std::string& SampleName, std::map< std::string, std::string >& Genome_Table, std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result, int barcode_index=0 )
	{
		this_analyzer_count_ = this_analyzer_count;
/*Joye
		if( this_analyzer_count == 0 )
		{
			/// @brief input, filter = each, db_idx = -1, db_depth=0, db_depth_value, sys
			Analysis (in, -1, -1, 0, "-2", 0);
			Analysis (in, 1, -1, 0, "-2", 0);
		}
*/
		Analysis( in, FilterType::value, DbIndexType::value, DbDepthType::value, boost::mpl::c_str<DbDepthNameType>::value, Genome_Table, 0 );//Joye
//Joye		Analysis( in, FilterType::value, DbIndexType::value, DbDepthType::value, boost::mpl::c_str<DbDepthNameType>::value );

		{
			std::lock_guard< std::mutex > lock( gOutMutex_ );

			if( gOutSet_.size() < this_analyzer_count + 1 )
			{
				gOutSet_.resize( this_analyzer_count + 1 );
				gOutNamePrefix_.resize( this_analyzer_count + 1 );
				//std::cerr << "this_analyzer_count " << this_analyzer_count << std::endl;
			}

			MapSum( gOutSet_[this_analyzer_count_], out_set_ );
			//std::cerr << "gOutSet_[this_analyzer_count_] *" << &gOutSet_[this_analyzer_count_] << std::endl;
		}
		//std::cerr << "this_analyzer_count_ " << this_analyzer_count_ << std::endl;

		//if 最後一個步驟
		if( eof_flag )
		{
			std::string db_depth_name_str( boost::mpl::c_str<DbDepthNameType>::value );
			std::string db_depth2_name_str( boost::mpl::c_str<DbDepth2NameType>::value );

			db_depth_name_str = convert_name( db_depth_name_str, 1 );
			db_depth2_name_str = convert_name( db_depth2_name_str, 2 );
 
			//lendist.Normal.full.GMPM
			std::string filename_prefix( "mirdist." );
			filename_prefix += boost::mpl::c_str< typename GetReadSeqClass::NAME_TYPE >::value;
			filename_prefix += ".";
			filename_prefix += boost::mpl::c_str< typename GetReadLengthClass::NAME_TYPE >::value;
			filename_prefix += ".";
			filename_prefix += boost::mpl::c_str< typename CalReadCountClass::NAME_TYPE >::value;
			
			filename_prefix += ".";
			filename_prefix += db_depth_name_str;
			filename_prefix += ".";
			filename_prefix += db_depth2_name_str;

			gOutNamePrefix_[this_analyzer_count] = filename_prefix; 	
		}

		return (void*) &in;//NULL;
	}

	std::string convert_name( std::string str, int type )
	{
		if( type == 1 )
		{
			if( str == "-1" )
				return "biotype";
			else if( str == "-2" )
				return "merged";
			else if( str == "-3" )
				return "no_biotype";
			else
				return str;
		}
		else
		{
			if( str == "-1" )
				return "detail";
			else if( str == "-2" )
				return "merged_detail";
			else if( str == "-3" )
				return "no_detail";
			else
				return str;
		}
	}

	//printer() has been delete
	
	/**
	 * @brief Analysis在run時，所呼叫的function，主要是判斷此 read 是不是要記錄
	 * @param is_filter Filter後 read(anno_rawbed)內部所記錄的 filter tag
	 * @param set_filter 使用者參數所設定的 filter tag (要filter 哪個 tag)，配合 is_filter 來決定 read 是否納入計算。
	 *		  -1 表示全部納入計算，0 表示去掉 tag=0，1表示去掉tag=1
	 * @return bool 回傳是否 filter掉（不納入計算）
	 */
	inline bool IsFilter( const int is_filter, const int set_filter )
	{
		if( set_filter == -1 )
			return false;
		else
		{
			if( is_filter == set_filter )
				return true;
			else
				return false;
		}
	}

	//Analysis2() has been delete	
	
	/**
	 * @brief Analysis實作，對 read做計算，記入進 output
	 * @param in 輸入的資料
	 * @param filter 使用者定的 filter 參數，詳情請看上面的 typedef
	 * @param db_index 使用者定的 db_index 參數，詳情請看上面的 typedef
	 * @param db_depth 使用者定的 db_depth 參數，詳情請看上面的 typedef
	 * @param db_depth_name 使用者定的 db_depth_name 參數，詳情請看上面的 typedef
	 * @return void
	 */
	void Analysis( INPUT_TYPE in, const int filter, const int db_index, const int db_depth, const char* db_depth_name, std::map< std::string, std::string >& Genome_Table, int sys = 2 )//Joye
//Joye	void Analysis( INPUT_TYPE in, const int filter, const int db_index, const int db_depth, const char* db_depth_name, int sys = 2 )
	{
		for( auto &anno_rawbed : *in )
		{
			if( IsFilter( anno_rawbed.is_filtered_, filter ))
				continue;
			
			for( int db_idx(0); db_idx != anno_rawbed.annotation_info_.size(); ++db_idx )
			{
				int check_db_idx = get_key( db_index, db_idx );
				if( check_db_idx == -3 )
					continue;

				for( int db_dep(0); db_dep != anno_rawbed.annotation_info_[db_idx].size(); ++db_dep )
				{
					int check_db_depth = get_key( db_depth, db_dep );
					if(check_db_depth == -3)
						continue;
					
					//設定的 db_depth_name 文字轉 string
					std::string db_depth_name_str( db_depth_name );

					//取得 annotation 內容
					std::string &db_depth_value = anno_rawbed.annotation_info_[db_idx][db_dep];

					//轉換 annotation 內容，db_depth_name_str 為 -1(全部分開), -2(全部合併), -3(不做), other(指定)
					std::string check_depth_name = get_key( db_depth_name_str, db_depth_value );

					if( check_depth_name == "-3" )
						continue;
					
					//設定的 db_depth2_name 文字轉 string
					std::string db_depth2_name_str( boost::mpl::c_str<DbDepth2NameType>::value );
					
					// for miRNA，如果db_depth2_name_str != -3 ，改變 check_depth_name，也就是改變記錄的文字，主要是可以讓第二層文字取代
					// e.g., 上面先判斷是不是 miRNA，然後如果 db_depth2_name_str = -1，那就會把 miRNA 的所有細分名字做 length distribution

					if( db_depth2_name_str != "-3" )
					{
						std::string &db_depth2_value = anno_rawbed.annotation_info_[db_idx][1];
						std::string check_depth2_name = get_key( db_depth2_name_str, db_depth2_value );
						check_depth_name = check_depth2_name;
					}	
					
					CalOutput< decltype( anno_rawbed )> 
					( anno_rawbed, sys, filter, check_db_idx, check_db_depth, check_depth_name.c_str(), Genome_Table );//Joye
//Joye					( anno_rawbed, sys, filter, check_db_idx, check_db_depth, check_depth_name.c_str() );
				}
			}
		}
	}
	
	std::string get_key( const std::string &type, const std::string &value )
	{
		if( type == "-1" )
		{
			return value;
		}
		else if( type == "-2" )
		{
			return "";
		}
		else
		{
			if( type == value )
			{
				return value;
			}
			else
			{
				return "-3";
			}
		}
		
	}
	int get_key( const int &type, const int &value )
	{
		if( type == -1 )
		{
			return value;
		}
		else if( type == -2 )
		{
			return -2;
		}
		else
		{
			if(type == value)
			{
				return value;
			}
			else
			{
				return -3;
			}
		}
	}
	
	
	/**
	 * @brief Analysis在run時，所呼叫的function，主要是簡化重複程式碼，此部份主要負責 insert into map of output
	 * @tparam READ_TYPE 自動決定 anno rawbed data type
	 * @param anno_rawbed read 的資料結構，詳情請看 annotation_raw_bed.hpp
	 * @param sys 此為特殊參數，主要是output所必須記錄，0=> db_index為空 and db_depth為空，1=> db_index為空 and db_depth不為空，2=>db_index不為空 and db_depth不為空
	 * @param filter 使用者定的 filter 參數，詳情請看上面的 typedef
	 * @param db_index 使用者定的 db_index 參數，詳情請看上面的 typedef
	 * @param db_depth 使用者定的 db_depth 參數，詳情請看上面的 typedef
	 * @param db_depth_name 使用者定的 db_depth_name 參數，詳情請看上面的 typedef
	 * @return void
	 */
	template<class READ_TYPE>
	inline void CalOutput( READ_TYPE &anno_rawbed ,const int sys, const int filter, const int db_index, const int db_depth, const char* db_depth_name, std::map< std::string, std::string >& Genome_Table)//Joye
//Joye	inline void CalOutput( READ_TYPE &anno_rawbed ,const int sys, const int filter, const int db_index, const int db_depth, const char* db_depth_name )
	{
		//std::cerr << "CalOutput: " << filter << "\t" << db_index << "\t" << db_depth << "\t" << db_depth_name << std::endl;
		//std::cerr << "----------------------------\nAnnoRawBed: " << anno_rawbed << "\n----------------------------\n";
		
		int read_length = GetReadLengthClass::GetReadLength( anno_rawbed );
		double count = CalReadCountClass::CalReadCount( anno_rawbed );
		std::string read_seq = GetReadSeqClass::GetReadSeq(anno_rawbed, Genome_Table);//Joye
//Joye		std::string read_seq = GetReadSeqClass::GetReadSeq( anno_rawbed );
		
//Joye		std::get<0>( out_set_[-1][std::make_tuple( sys,filter,db_index,db_depth,db_depth_name, read_seq )]) += count;
		std::get<0>( out_set_[-1][std::make_tuple( sys,filter,db_index,db_depth,"", "" )]) += count; //Joye
		std::get<0>( out_set_[read_length][std::make_tuple( sys,filter,db_index,db_depth,db_depth_name, read_seq )]) += count;
	}
	
	void MapSum( OutPutType &map_sum, OutPutType &map_item )
	{
		map_sum = map_sum + map_item;
		/*
		//第一層 map<len, map<tuple, tuple> >
		for(auto &map_item_pair : map_item)
		{
			//第二層 map<tuple,tuple>
			auto &map_item_pair_item =  map_item_pair.second;
			
			for(auto &map_item_pair_item_pair : map_item_pair_item)
			{
				std::get<0>(map_sum[map_item_pair.first][map_item_pair_item_pair.first]) += std::get<0>(map_item_pair_item_pair.second);
			}
			
		}
		*/
	}
	
	
};


template<class INPUT_TYPE>
std::vector< std::map<int, std::map<std::tuple<int,int,int,int,std::string, std::string>, std::tuple<double, double, double> > > > 
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutSet_(0);

template<class INPUT_TYPE>
std::vector< std::string > 
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutNamePrefix_(0);

template<class INPUT_TYPE>
std::mutex
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutMutex_;


template<class INPUT_TYPE, class ANALYZER_TYPELIST>
std::vector< std::map<int, std::map<std::tuple<int,int,int,int,std::string, std::string>, std::tuple<double, double, double> > > > &
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>::gOutSet_
= AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutSet_;

template<class INPUT_TYPE, class ANALYZER_TYPELIST>
std::vector< std::string > &
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>::gOutNamePrefix_
= AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutNamePrefix_;

template<class INPUT_TYPE, class ANALYZER_TYPELIST>
std::mutex &
AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerTypes::miRNADistribution>::gOutMutex_
= AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST_GLOBAL_MIRDIST, AnalyzerTypes::miRNADistribution>::gOutMutex_;

#endif
