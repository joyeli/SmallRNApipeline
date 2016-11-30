#ifndef ANALYZER_POLICY_HPP_
#define ANALYZER_POLICY_HPP_

#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <type_traits>
#include <boost/type_traits.hpp>
#include <boost/mpl/char.hpp>
#include <boost/mpl/string.hpp>

#include <pokemon/constant_def.hpp>
#include <pokemon/tuple_utility.hpp>
#include <pokemon/analyzer/analyzer_utility.hpp>

struct GetReadFullLength
{
	typedef boost::mpl::string<'full'> NAME_TYPE;
	
	template<class READ_TYPE>
	static int GetReadLength(READ_TYPE &raw_read)
	{
		return raw_read.length_;
	}
};
struct GetReadPrefixLength
{
	typedef boost::mpl::string<'pre', 'fix'> NAME_TYPE;
	
	template<class READ_TYPE>
	static int GetReadLength(READ_TYPE &raw_read)
	{
		return raw_read.length_ - raw_read.tail_length_;
	}
};
struct CalReadCountGMPM
{
	typedef boost::mpl::string<'GMPM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		return (double) raw_read.reads_count_ / raw_read.multiple_alignment_site_count_;
	}
};
struct CalReadCountGMOnly
{
	typedef boost::mpl::string<'GM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		if(raw_read.tail_length_ == 0)
			return (double) raw_read.reads_count_ / raw_read.multiple_alignment_site_count_;
		else
			return 0.0;
	}
};
struct CalReadCountPMOnly
{
	typedef boost::mpl::string<'PM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		if(raw_read.tail_length_ > 0)
			return (double) raw_read.reads_count_ / raw_read.multiple_alignment_site_count_;
		else
			return 0.0;
	}
};
struct CalReadSpeciesGMPM
{
	typedef boost::mpl::string<'GMPM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		return ((double) 1 / raw_read.multiple_alignment_site_count_);
	}
};
struct CalReadSpeciesGMOnly
{
	typedef boost::mpl::string<'GM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		if(raw_read.tail_length_ == 0)
			return (double) 1 / raw_read.multiple_alignment_site_count_;
		else
			return 0.0;
	}
};
struct CalReadSpeciesPMOnly
{
	typedef boost::mpl::string<'PM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		if(raw_read.tail_length_ > 0)
			return (double) 1 / raw_read.multiple_alignment_site_count_;
		else
			return 0.0;
	}
};




/**
 * @brief 回傳 read length。取得 anno raw bed length 的方法，也就是length distribution 的 length
 */
struct GetReadLengthDefault
{
	typedef boost::mpl::string<'full'> NAME_TYPE;
	
	template<class READ_TYPE>
	static int GetReadLength(READ_TYPE &raw_read)
	{
		return raw_read.length_;
	}
};

//============Policy for CalReadCount=================

/**
 * @brief 回傳 read count。取得 anno raw bed count 的方法，也就是length distribution 的 value
 */
struct CalReadCountDefault
{
	typedef boost::mpl::string<'GMPM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		//std::cerr<<"GetReadSeq of CalReadCount version return "<<(double)(raw_read.reads_count_ / raw_read.multiple_alignment_site_count_)<<std::endl;
		return (double) raw_read.reads_count_ / raw_read.multiple_alignment_site_count_;
	}
};

/**
 * @brief 回傳 species count。取得 anno raw bed count 的方法，也就是call species應用的 value
 */
struct CalReadSpecies
{
	typedef boost::mpl::string<'GMPM'> NAME_TYPE;
	
	template<class READ_TYPE>
	static double CalReadCount(READ_TYPE &raw_read)
	{
		//std::cerr<<"GetReadSeq of CalRadSpecies version return "<<(double)(1 / raw_read.multiple_alignment_site_count_)<<std::endl;
		return ((double) 1 / raw_read.multiple_alignment_site_count_);
	}
};

//====================================================

//============Policy for CalReadSeq===================

/**
 * @brief 回傳 read sequence。取得 anno raw bed sequence content 的方法
 */
struct GetReadSeqDefault
{
	typedef boost::mpl::string<'Norm', 'al'> NAME_TYPE;
	
	template<class READ_TYPE>
	static std::string GetReadSeq(READ_TYPE &raw_read)
	{
		return "";
	}
};

/**
 * @brief 回傳 read sequence。取得 anno raw bed tailing content 的方法
 */
struct GetReadTailing
{
	typedef boost::mpl::string<'Tail'> NAME_TYPE;
	
	template<class READ_TYPE>
	static std::string GetReadSeq(READ_TYPE &raw_read)
	{
		std::map <int, char> mtable ({ {0,'A'}, {1,'C'}, {2,'G'}, {3,'T'} });
		std::string temp;
		for (auto index=raw_read.tail_length_-1; index>=0; --index)
			temp.push_back ( mtable [(raw_read.tail_mismatch_ >> (2*index) ) % 4] );
		//std::cerr<<"raw_read_tail_mismatch_ & seq "<<raw_read.tail_mismatch_<<'\t'<<temp<<'\n';
		return temp;
	}
};

/**
 * @brief 回傳 read sequence。取得 anno raw bed sequence content 的方法
 */
template <int FirstOrLast, int START, int LENGTH>
struct GetReadFirstNLastComposition
{
	/// @brief First-1nt, First-2nt, Last-1nt, Last-2nt
	typedef typename boost::mpl::if_
	<
		typename boost::is_same< boost::mpl::int_<FirstOrLast>, boost::mpl::int_<0> >
		,boost::mpl::string<'Firs', 't-', LENGTH+48, 'nt'>
		,boost::mpl::string<'Las', 't-', LENGTH+48, 'nt'>
	>::type NAME_TYPE;
	
	
	template<class READ_TYPE>
	static std::string GetReadSeq(READ_TYPE &raw_read, std::map< std::string, std::string >& Genome_Table)//Joye
//Joye	static std::string GetReadSeq(READ_TYPE &raw_read)
	{
		/// @breif First = 0
		if ( ! FirstOrLast )
			return GetReadContent (raw_read, START, LENGTH, Genome_Table);//Joye
//Joye			return GetReadContent (raw_read, START, LENGTH);
		/// @breif Last != 0
		else
			return GetReadContent (raw_read, raw_read.length_ - START - LENGTH, LENGTH, Genome_Table);//Joye
//Joye			return GetReadContent (raw_read, raw_read.length_ - START - LENGTH, LENGTH);
	}

	template<class READ_TYPE>
	static std::string GetReadContent (READ_TYPE &target, const int start, const int length, std::map< std::string, std::string >& Genome_Table)//Joye
//Joye	static std::string GetReadContent (READ_TYPE &target, const int start, const int length)
	{
/*
		std::stringstream chrstream;
		chrstream <<target.chr_prefix_;

        if( target.chr_idx_>=1 && target.chr_idx_<=32 )
            chrstream << (int)target.chr_idx_;
        else if( target.chr_idx_>=33 && target.chr_idx_<=64 )
            chrstream << (int)target.chr_idx_-32;
        else if( target.chr_idx_>=65 && target.chr_idx_<=90 )
            chrstream << target.chr_idx_;
        else if( target.chr_idx_>=97 && target.chr_idx_<=122)
            chrstream << (char)(target.chr_idx_ -32);

		std::string chr = chrstream.str();

*/
		//std::string chr = target.GetChr (target.chr_idx_);
		//std::string return_string = (PipelinePreparator<>::gGenometable_[chr]).substr (target.start_+start, length);
		//std::transform(return_string.begin(), return_string.end(),return_string.begin(), ::toupper );
		std::string return_string = target.getReadSeq(Genome_Table);//Joye
//Joye		std::string return_string = target.getReadSeq();
		
		return return_string.substr(start, length);
		
		//return std::string (PipelinePreparator<>::gGenometable_[chr]).substr (target.start_+start, length);
	}
	
};

/**
 * @brief 回傳 read sequence。取得 anno raw bed tailing content 的方法
 */
template < int START=1, int LENGTH=6 >
struct GetReadSeed
{
	typedef boost::mpl::string<'seed'> NAME_TYPE;
	
	template<class READ_TYPE>
	static std::string GetReadSeq(READ_TYPE &raw_read, std::map< std::string, std::string >& Genome_Table)//Joye
//Joye	static std::string GetReadSeq(READ_TYPE &raw_read,)
	{
		return GetReadContent (raw_read, START, LENGTH, Genome_Table);//Joye
//Joye		return GetReadContent (raw_read, START, LENGTH);
	}

	template<class READ_TYPE>
	static std::string GetReadContent (READ_TYPE &target, const int start, const int length, std::map< std::string, std::string >& Genome_Table)//Joye
//Joye	static std::string GetReadContent (READ_TYPE &target, const int start, const int length)
	{
/*		std::stringstream chrstream;
		chrstream <<target.chr_prefix_;

        if( target.chr_idx_>=1 && target.chr_idx_<=32 )
            chrstream << (int)target.chr_idx_;
        else if( target.chr_idx_>=33 && target.chr_idx_<=64 )
            chrstream << (int)target.chr_idx_-32;
        else if( target.chr_idx_>=65 && target.chr_idx_<=90 )
            chrstream << target.chr_idx_;
        else if( target.chr_idx_>=97 && target.chr_idx_<=122)
            chrstream << (char)(target.chr_idx_ -32);

		std::string chr = chrstream.str();

*/
		//std::string chr = target.GetChr (target.chr_idx_);
		//return std::string//("");
		//		(PipelinePreparator<>::gGenometable_[chr]).substr (target.start_+start, length);
		std::string return_string = target.getReadSeq(Genome_Table);//Joye
//Joye		std::string return_string = target.getReadSeq();
		return return_string.substr(start, length);
	}
};

//====================================================




#endif
