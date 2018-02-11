/**
 *  @file analyzer.hpp
 *  @brief analyzer 計算 filter 後的結果
 *  @author C-Salt Corp.
 */
#ifndef ANALYZER_HPP_
#define ANALYZER_HPP_

#include <string>
#include <iostream>
#include <vector>
#include <tuple>
#include <mutex>
#include <set>

#include <boost/type_traits.hpp>
#include <boost/type_traits/add_pointer.hpp>
#include <boost/mpl/char.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/for_each.hpp>



#include <pokemon/constant_def.hpp>
#include <pokemon/tuple_utility.hpp>
//#include <pokemon/pipeline/pipeline_preparator.hpp>

#include <pokemon/analyzer/analyzer_utility.hpp>
#include <AGO/algorithm/analyzer_policy_ForAgoSorting.hpp>


/**
 * @brief Analyzer 的實作，此為範型版本，用來特化
 * @tparam INPUT_TYPE 輸入資料的型別，一定為 vector，通常為 vector<AnnoRawBed>
 * @tparam ANALYZER_TYPELIST Analyzer要用的參數設定，通常為 boost::mpl::vector<boost::mpl::map<boost::mpl::pair<KEY, VALUE> > >
 * @tparam ANALYZER_TYPE 特化 Analyzer用的參數，為constant_def 內的 AnalyzerTypes
 */
template<class INPUT_TYPE, class ANALYZER_TYPELIST, int ANALYZER_TYPE>
class AnalyzerImpl
{};

/**
 * @brief Analyzer是以type list的形態出現，此class為Analyzer boost::mpl::for_each所直接執行的，可以當做一個中介產物
 * @tparam INPUT_TYPE 通常為vector<AnnoRawBed>
 */
template<class INPUT_TYPE>
class AnalyzerImplInitC
{
public:
	INPUT_TYPE in;
	
	/// @brief pipe_index_ 第N輪呼叫的次序，原始數據會被分群後的順序。目的是要在平行跑完後，最終能夠還原順序，是由外面傳進來的
    size_t pipe_index_;
    
    /// @brief eof_flag_ 是否為最後一輪（也就是最後一次呼叫），由外面傳進來 
    bool eof_flag_;	
   
	std::string& OutputPath_; //Joye
	std::string& SampleName_; //Joye
	std::map< std::string, std::string >& Genome_Table_; //Joye
	std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result_; //Joye
 
    /// @brief barcode_index_ smaple在分完Barcode後，所記錄的barcode次序
    size_t barcode_index_;
    
    //keep track of information of how many times, the mapped value, a certain kind of AnalyzerType, the key value, has been involved in the current ANALYZER_TYPELIST
	/// @brief analyzer_count_type_ ，key為不同的 analyzer，value為某analyzer被run了幾次（這邊與第幾輪無關，僅僅是pipeline的順序次數）
	std::map<int,int> analyzer_count_type_; 
	
	/// @brief Global儲存的指標，主要是可以讓外面取得到特定物件
	static void* gPtr_;
	
	AnalyzerImplInitC(INPUT_TYPE i)
		: in (i)
	{}

	~AnalyzerImplInitC()
	{}

	AnalyzerImplInitC()
	{}

	AnalyzerImplInitC(INPUT_TYPE i, size_t pipe_index, bool eof_flag, std::string& OutputPath, std::string& SampleName, std::map< std::string, std::string >& Genome_Table, std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result, size_t barcode_index) //Joye
        : in (i), pipe_index_ (pipe_index), eof_flag_ (eof_flag), OutputPath_ (OutputPath), SampleName_ (SampleName), Genome_Table_ (Genome_Table), Analyzer_Result_ (Analyzer_Result), barcode_index_ (barcode_index) //Joye
//Joye	AnalyzerImplInitC(INPUT_TYPE i, size_t pipe_index, bool eof_flag, std::string OutputPath, std::string SampleName, size_t barcode_index)
//Joye        : in (i), pipe_index_ (pipe_index), eof_flag_ (eof_flag), barcode_index_ (barcode_index)
    {}

	/// @tparam ANALYZER_TYPELIST type list 通常型別為 boost::mpl::map，負責記錄使用者訂定的 analyzer
	template <class ANALYZER_TYPELIST>
	void operator()(ANALYZER_TYPELIST t)
	{
		/// @breif 取得 analyzer type (enum)
		typedef typename boost::mpl::at<ANALYZER_TYPELIST, boost::mpl::int_<0> >::type AnalyzerType;
		/// @brief 建構 analyzer
		AnalyzerImpl<INPUT_TYPE, ANALYZER_TYPELIST, AnalyzerType::value> Analysis(in);
//		gPtr_ = 
		Analysis(analyzer_count_type_[AnalyzerType::value], pipe_index_, eof_flag_, OutputPath_, SampleName_, Genome_Table_, Analyzer_Result_, barcode_index_);
//Joye		Analysis(analyzer_count_type_[AnalyzerType::value], pipe_index_, eof_flag_, barcode_index_);	
		++analyzer_count_type_[AnalyzerType::value];
	}
};

/// @brief 初始化 static 物件
template<class INPUT_TYPE>
void*
AnalyzerImplInitC<INPUT_TYPE>::gPtr_ = (void*) NULL;


/**
 * @brief Analyzer 的主要介面，要將 type list 做 boost mpl for_each，依序把不同的 analyzer做完
 * @tparam INPUT_TYPE 通常為vector<AnnoRawBed>
 * @tparam ANALYZER_TYPELIST type list 通常型別為 boost::mpl::map，負責記錄使用者訂定的 analyzer
 */
template<class INPUT_TYPE, class ANALYZER_TYPELIST>
class AnalyzerC
{
public:
	~AnalyzerC()
	{}
	
	/**
	 * @brief 外面可以呼叫此介面，開始執行 analyzer
	 * @param pipe_index_ 第N輪呼叫的次序，原始數據會被分群後的順序。目的是要在平行跑完後，最終能夠還原順序
	 * @param eof_flag_ 是否為最後一輪（也就是最後一次呼叫）
	 * @return 回傳 AnalyzerImplInit 的static指標，主要是給予外面使用一些物件的能力
	 */
	void* run (INPUT_TYPE in, size_t pipe_index, bool eof_flag) 
	{
		boost::mpl::for_each<ANALYZER_TYPELIST> ( AnalyzerImplInitC<INPUT_TYPE>( in, pipe_index, eof_flag) );
		return AnalyzerImplInitC<INPUT_TYPE>::gPtr_;
	}

	void* run (INPUT_TYPE in, size_t pipe_index, bool eof_flag, std::string& OutputPath, std::string& SampleName, std::map< std::string, std::string >& Genome_Table, std::vector< std::map< std::string, std::map< std::string, double >>>& Analyzer_Result, size_t barcode_index ) //Joye
//Joye	void* run (INPUT_TYPE in, size_t pipe_index, bool eof_flag, size_t barcode_index ) 
	{
		boost::filesystem::path dir( OutputPath );
//Joye		boost::filesystem::path dir("output/");
		boost::filesystem::create_directory(dir);
//		auto output_path = std::string("output/sample-") + std::to_string(barcode_index) + std::string("/");
//Joye		auto output_path = std::string("output/sample-") + PipelinePreparator<>::gBarcode_vector_[barcode_index] + std::string("/");
		std::string output_path = OutputPath + "/" + SampleName + "/"; //Joye
		boost::filesystem::path dir_analyzer(output_path);
		boost::filesystem::create_directory(dir_analyzer); 

		boost::mpl::for_each<ANALYZER_TYPELIST> ( AnalyzerImplInitC<INPUT_TYPE>( in, pipe_index, eof_flag, OutputPath, SampleName, Genome_Table, Analyzer_Result, barcode_index) ); //Joye
//Joye		boost::mpl::for_each<ANALYZER_TYPELIST> ( AnalyzerImplInitC<INPUT_TYPE>( in, pipe_index, eof_flag, barcode_index) );
		return AnalyzerImplInitC<INPUT_TYPE>::gPtr_;
	}
};

/**
 * @brief Analyzer 有許多不同的功能，此 class 是定義不同功能，所提供的不同參數。目的是讓使用者在使用參數時，要指定一個參數名，方便理解程式
 * @tparam ANALYZER_TYPE constant_def 內 AnalyzerTypes 的 type，用來特化不同的 Analyzer功能
 */

template<int ANALYZER_TYPE>
class AnalyzerParameter
{
public:
	/** 
	 * @brief 定義參數有哪些
	 */
	/// @brief AnalyzerType 特化Analyzer的type，每一個AnalyzerParameter都應該要有此參數
	typedef boost::mpl::int_<0> AnalyzerType;
};


#include "analyzer_length_distribution_ForAgoSorting.hpp"
#include "analyzer_length_distribution_printer_ForAgoSorting.hpp"
#include "analyzer_mirna_distribution_ForAgoSorting.hpp"
#include "analyzer_mirna_distribution_printer_ForAgoSorting.hpp"

#endif
