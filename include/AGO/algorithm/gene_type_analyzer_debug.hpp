#pragma once
#include <AGO/algorithm/gene_type_analyzer_declare.hpp>

namespace ago {
namespace algorithm {

void debug_outputing( std::vector< std::vector< CountingTableType >>& anno_table_tail )
{
    for( std::size_t smp = 0; smp < anno_table_tail.size(); ++smp )
    {
        for( std::size_t tail = 0; tail < anno_table_tail[ smp ].size(); ++tail )
        {
            for( auto& anno : anno_table_tail[ smp ][ tail ])
            {
                for( auto& len : anno.second )
                    std::cerr << std::setprecision( 2 ) << std::fixed
                        << smp << "\t"
                        << tail << "\t"
                        << anno.first << "\t"
                        << len.first << "\t"
                        << len.second << "\n"; 
            }
        }
    }
};

} // end of namespace algorithm
} // end of namespace ago
