#pragma once
#include <fstream>
#include <algorithm>
#include <string>
#include <cmath>
#include <tuple>
#include <map>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <vector>
#define Vector std::vector

template< class ValueType = Vector< double >>
class TmmDataType
{
public :

    ValueType values_;
    ValueType Mgrk_;
    ValueType Wgrk_;
    ValueType Mg_;
    ValueType Ag_;
    ValueType Yg_;
    double Nk_;
    double tmm_;

    TmmDataType()
        : values_( ValueType() )
        , Mgrk_( ValueType() )
        , Wgrk_( ValueType() )
        , Mg_( ValueType() )
        , Ag_( ValueType() )
        , Nk_( double() )
        , tmm_( double() )
    {}

    TmmDataType( ValueType& value )
        : values_( ValueType() )
        , Mgrk_( ValueType() )
        , Wgrk_( ValueType() )
        , Mg_( ValueType() )
        , Ag_( ValueType() )
        , Nk_( double() )
        , tmm_( double() )
    {
        values_ = std::move( value );
        Nk_ = NkCalculation();
    }

    double NkCalculation()
    {
        double Nk = 0;
        for( auto& value : values_ )
            Nk = Nk + value;
        return Nk;
    }

    void MAgCalculation()
    {
        double value_div_nk = 0.0;
        double trim_mean_div_nk = 1.0 / double( values_.size() );

        for( auto& value : values_ )
        {
            value_div_nk = value / Nk_;
            Mg_.push_back(  log2(( value_div_nk )/( trim_mean_div_nk )));
            Ag_.push_back(( log2(( value_div_nk )*( trim_mean_div_nk )))/2 );
        }
    }

    void ValueTrimming( auto& trim_mean, auto& MgPer, auto& AgPer )
    {
        double trim_nk = trim_mean * values_.size();
        double log2_trim_mean_div_nk = log2( 1.0 / double( values_.size() ));
        double trim_temp = (( trim_nk - trim_mean )/( trim_nk * trim_mean ));

        double mg_rate = double(Mg_.size()) / 100;
        double mg_lowr = MgPer * mg_rate;
        double mg_uppr = Mg_.size() - mg_lowr;

        double ag_rate = double(Ag_.size()) / 100;
        double ag_lowr = AgPer * ag_rate;
        double ag_uppr = Ag_.size() - ag_lowr;

        Vector< std::pair< double, size_t >> Mg, Ag;
        std::set< size_t > Mg_trim, Ag_trim;
        Vector< size_t > trimmed_vector;

        for( int i = 0; i < values_.size(); ++i )
        {
            Mg.push_back( std::make_pair( Mg_[i], i ));
            Ag.push_back( std::make_pair( Ag_[i], i ));
        }

        std::sort( Mg.begin(), Mg.end(),
            []( const std::pair< double, size_t >& a, const std::pair< double, size_t >& b )
            { return a.first < b.first; });

        std::sort( Ag.begin(), Ag.end(),
            []( const std::pair< double, size_t >& a, const std::pair< double, size_t >& b )
            { return a.first < b.first; });

        for( int i = 0; i < values_.size(); ++i )
        {
            if( mg_lowr <= i && mg_uppr >= i ) Mg_trim.emplace( Mg[i].second );
            if( ag_lowr <= i && ag_uppr >= i ) Ag_trim.emplace( Mg[i].second );
        }

        for( auto& Mi : Mg_trim )
            if( Ag_trim.find( Mi ) != Ag_trim.end() )
                trimmed_vector.push_back( Mi );

        for( auto& i : trimmed_vector )
        {
            if( values_[i] <= 0 ) continue;
            Mgrk_.push_back( log2( values_[i] / Nk_ ) / log2_trim_mean_div_nk );
            Wgrk_.push_back( (( Nk_ - values_[i] ) / ( Nk_ * values_[i] )) + trim_temp );
        }
    }

    void TmmCalculation()
    {
        double WMrk = 0;
        double Wrk  = 0;

        for( int i = 0; i < Mgrk_.size(); ++i )
        {
            WMrk = WMrk + ( Wgrk_[i] * Mgrk_[i] );
            Wrk = Wrk + Wgrk_[i];
        }

        tmm_ = WMrk / Wrk;
    }

    void TmmNormalization()
    {
        for( auto& value : values_ )
            value = value / tmm_;
    }
};


class TmmNor
{
public:

    TmmNor(
            Vector< TmmDataType<> >& samples,
            Vector< double >& tmm_vector,
            const double& Mg_per,
            const double& Ag_per,
            const std::size_t& trim_mean
    )
    {
        for( auto& sample : samples )
        {
            sample.MAgCalculation();
            sample.ValueTrimming( trim_mean, Mg_per, Ag_per );
            sample.TmmCalculation();
            sample.TmmNormalization();
            tmm_vector.emplace_back( sample.tmm_ );
        }
    }
};
