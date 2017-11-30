#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

template< class ValueType = std::vector< std::pair< double, int >>, class DrewType = std::vector< bool >>
class QuantileDataType
{
public :

	ValueType value_;
	DrewType  drew_;

	QuantileDataType()
	: value_( ValueType() )
	, drew_ ( DrewType() )
	{}

	QuantileDataType( std::vector< double >& value )
	: value_( ValueType() )
	, drew_ ( DrewType() )
	{
		for( int idx = 0; idx < value.size(); ++idx )
		{
			value_.emplace_back( std::make_pair( std::move(value[idx]), idx ));
			drew_.emplace_back( idx != 0 && value[idx] == value[idx-1] ? true : false );
		}
	}
};

class QuantileNor
{
public:
	QuantileNor( std::vector< QuantileDataType<> >& qdata )
	{
		//PrintQuData( qdata, "QuntileData" );

		RankSorting( qdata );
		//PrintQuData( qdata, "RankSorting" );

		RankMeanCal( qdata );
		//PrintQuData( qdata, "RankMeanCal" );

		ResumeOrder( qdata );
		//PrintQuData( qdata, "ResumeOrder" );

		ReplaceDrew( qdata );
		//PrintQuData( qdata, "ReplaceDrew" );
	}

	void RankSorting( auto& qdata )
	{
		for( auto& sample : qdata )
		{
			std::sort( sample.value_.begin(), sample.value_.end(),
					[]( const std::pair< double, int >& a, const std::pair< double, int >& b )
					{
						if( a.first == b.first )
                        {
							return a.second < b.second;
                        }
						else
                        {
							return a.first < b.first;
                        }
					});
		}
	}

	void RankMeanCal( auto& qdata )
	{
        double sum = 0.0;
		double mean = 0.0;

		for( int idx = 0; idx < qdata[0].value_.size(); ++idx )
		{
            sum = 0.0;
			for( auto& sample : qdata )
            {
				sum = sum + sample.value_[ idx ].first;
            }

            mean = sum / (double)(qdata.size());
			for( auto& sample : qdata )
            {
				sample.value_[ idx ].first = mean;
            }
		}
	}

	void ResumeOrder( auto& qdata )
	{
		for( auto& sample : qdata )
		{
			std::sort( sample.value_.begin(), sample.value_.end(),
					[]( const std::pair< double, int >& a, const std::pair< double, int >& b )
					{
						if( a.second == b.second )
                        {
							return a.first < b.first;
                        }
						else
                        {
							return a.second < b.second;
                        }
					});
		}
	}

	void ReplaceDrew( auto& qdata )
	{
		for( int idx = 0; idx < qdata[0].value_.size(); ++idx )
        {
			for( auto& sample : qdata )
            {
				if( sample.drew_[ idx ] == true )
                {
					sample.value_[ idx ].first = sample.value_[idx-1].first; 
                }
            }
        }
	}

	void PrintQuData( auto& qdata, auto stage )
	{
		std::cerr << "\n" << stage << ":\n";
		for( int i = 0; i < qdata[0].value_.size(); i++ )
		{
			for( auto& q : qdata )
            {
				std::cerr << "\t" << q.value_[i].first << "/" << q.value_[i].second << "/" << q.drew_[i];
            }

			std::cerr << "\n";
		}
	}
};
