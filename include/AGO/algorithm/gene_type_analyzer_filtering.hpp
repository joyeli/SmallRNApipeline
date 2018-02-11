#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerFiltering
{
  public:

    GeneTypeAnalyzerFiltering()
    {}

    using DropTypeList = boost::mpl::vector<
        boost::mpl::vector< boost::mpl::string< 'rm', 'sk' >, boost::mpl::int_< 1 >, boost::mpl::char_< '=' >>>;

    using Filters = FilterWorker< AnnotationRawBed<>, DropTypeList >;

    void drop_filtering( std::vector< BedSampleType >& bed_samples )
    {
		Filters run_filter;
        ParaThreadPool smp_parallel_pool( bed_samples.size() );
        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            smp_parallel_pool.job_post( [ smp, &run_filter, &bed_samples ] () mutable
            {
                for( auto& anno_rawbed : bed_samples[ smp ].second )
                {
                    anno_rawbed = run_filter.Filter( anno_rawbed );
                }
            });
        }
        smp_parallel_pool.flush_pool();
    }

};

} // end of namespace algorithm
} // end of namespace ago
