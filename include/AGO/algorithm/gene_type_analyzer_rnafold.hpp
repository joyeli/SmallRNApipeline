#pragma once
#include <cstdlib>

namespace ago {
namespace algorithm {

class GeneTypeAnalyzerRNAfold
{

  public:

    GeneTypeAnalyzerRNAfold()
    {}

    static std::vector< std::map< std::string, std::tuple< double, double, double, double >>> make_fold_table(
            std::vector< BedSampleType >& bed_samples,
            AnnoLengthIndexType& ano_len_idx,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::string& rnafold_path,
            const std::string& biotype,
            const std::string& token = ""
            )
    {
        if(( biotype != "miRNA_mirtron" && biotype != "miRNA" && biotype != "mirtron" ) || rnafold_path == "" )
            return std::vector< std::map< std::string, std::tuple< double, double, double, double >>>();
        
        std::vector< std::map< std::string, std::tuple< double, double, double, double >>> results( bed_samples.size() );
        //                                              5p      mid     3p      3pTailOnly
        
        for( auto& anno : ano_len_idx.first )
        {
            if( results[0].find( anno ) == results[0].end() )
            {
                for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
                    results[ smp ][ anno ] = { 0.0, 0.0, 0.0, 0.0 };
            }
        }
        
        std::set< std::string > annos;
        std::map< std::string, std::vector< std::string >> rnafolds;
        std::map< std::string, std::vector< std::map< std::string, double >>> entropies;
        std::map< std::string, std::vector< std::map< std::string, double >>> counts;

        GeneTypeAnalyzerSqalign::get_rnafold( rnafold_path, rnafolds );

        entropies[ "5p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "mid" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "3p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        entropies[ "3pTailOnly" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );

        counts[ "5p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "mid" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "3p"  ] = std::vector< std::map< std::string, double >>( bed_samples.size() );
        counts[ "3pTailOnly" ] = std::vector< std::map< std::string, double >>( bed_samples.size() );

        bool is_lens = token != "" ? true : false;

        std::string mir;

        std::string gene_name;
        std::string gene_seed;

        std::size_t strpos;
        std::size_t midpos;
        std::size_t endpos;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            for( auto& raw_bed : bed_samples[ smp ].second )
            {
                if( raw_bed.ppm_ < filter_ppm ) continue;
                if( raw_bed.annotation_info_.empty() || raw_bed.annotation_info_[0].empty() ) continue;
                for( int i = 0; i < raw_bed.annotation_info_[0].size(); i+=2 )
                {
                    if( raw_bed.annotation_info_[0][i] != "miRNA" &&
                        raw_bed.annotation_info_[0][i] != "mirtron" &&
                        raw_bed.annotation_info_[0][i] != "miRNA_mirtron" )
                        continue;

                    if( is_lens && token != std::to_string( raw_bed.length_ - raw_bed.tail_length_ ))
                        continue;

                    gene_name = raw_bed.annotation_info_[0][ i+1 ];
                    mir = gene_name.substr( 0, raw_bed.annotation_info_[0][ i+1 ].length() -3 );

                    if( raw_bed.start_ < std::stoll( rnafolds[ mir ][1] ) + 1 ) continue;
                    if( raw_bed.end_   > std::stoll( rnafolds[ mir ][2] ) + 1 ) continue;

                    if( entropies[ "5p"  ][ smp ].find( gene_name ) == entropies[ "5p"  ][ smp ].end() ) entropies[ "5p"  ][ smp ][ gene_name ] = 0.0;
                    if( entropies[ "mid" ][ smp ].find( gene_name ) == entropies[ "mid" ][ smp ].end() ) entropies[ "mid" ][ smp ][ gene_name ] = 0.0;
                    if( entropies[ "3p"  ][ smp ].find( gene_name ) == entropies[ "3p"  ][ smp ].end() ) entropies[ "3p"  ][ smp ][ gene_name ] = 0.0;

                    if( counts[ "5p"  ][ smp ].find( gene_name ) == counts[ "5p"  ][ smp ].end() ) counts[ "5p"  ][ smp ][ gene_name ] = 0.0;
                    if( counts[ "mid" ][ smp ].find( gene_name ) == counts[ "mid" ][ smp ].end() ) counts[ "mid" ][ smp ][ gene_name ] = 0.0;
                    if( counts[ "3p"  ][ smp ].find( gene_name ) == counts[ "3p"  ][ smp ].end() ) counts[ "3p"  ][ smp ][ gene_name ] = 0.0;

                    switch( raw_bed.strand_ )
                    {
                        case '+' : strpos = raw_bed.start_ - ( std::stoll( rnafolds[ mir ][1] ) + 1 ); break;
                        case '-' : strpos = ( std::stoll( rnafolds[ mir ][2] ) + 1 ) - raw_bed.end_  ; break;
                    }

                    annos.emplace( gene_name );
                    endpos = strpos + 1 + (int)raw_bed.length_ - (int)raw_bed.tail_length_;
                    std::vector< std::string > entropy( rnafolds[ mir ].begin() + 6, rnafolds[ mir ].begin() + rnafolds[ mir ].size() );

                    entropies[ "5p"  ][ smp ][ gene_name ] += get_entropy( strpos,     8, entropy );
                    entropies[ "mid" ][ smp ][ gene_name ] += get_entropy( strpos + 8, 4, entropy );
                    entropies[ "3p"  ][ smp ][ gene_name ] += get_entropy( endpos - 8, 8, entropy );

                    counts[ "5p"  ][ smp ][ gene_name ] += 1;
                    counts[ "mid" ][ smp ][ gene_name ] += 1;
                    counts[ "3p"  ][ smp ][ gene_name ] += 1;

                    if( GeneTypeAnalyzerCounting::seqT2U( raw_bed.getTail() ) != "" )
                    {
                        if( entropies[ "3pTailOnly" ][ smp ].find( gene_name ) == entropies[ "3pTailOnly" ][ smp ].end() ) entropies[ "3pTailOnly" ][ smp ][ gene_name ] = 0.0;
                        if( counts[ "3pTailOnly" ][ smp ].find( gene_name ) == counts[ "3pTailOnly" ][ smp ].end() ) counts[ "3pTailOnly" ][ smp ][ gene_name ] = 0.0;

                        entropies[ "3pTailOnly"  ][ smp ][ gene_name ] += get_entropy( endpos - 8, 8, entropy );
                        counts[ "3pTailOnly"  ][ smp ][ gene_name ] += 1;
                    }
                }
            }
        }

        for( auto& type : entropies )
        {
            for( std::size_t smp = 0; smp < type.second.size(); ++smp )
            {
                for( auto& anno : type.second[ smp ] )
                {
                    if( type.first == "5p"          ) std::get<0>( results[ smp ][ anno.first ] ) = anno.second / counts[ type.first ][ smp ][ anno.first ];
                    if( type.first == "mid"         ) std::get<1>( results[ smp ][ anno.first ] ) = anno.second / counts[ type.first ][ smp ][ anno.first ];
                    if( type.first == "3p"          ) std::get<2>( results[ smp ][ anno.first ] ) = anno.second / counts[ type.first ][ smp ][ anno.first ];
                    if( type.first == "3pTailOnly"  ) std::get<3>( results[ smp ][ anno.first ] ) = anno.second / counts[ type.first ][ smp ][ anno.first ];
                }
            }
        }

        return results;
    }

    static double get_entropy(
            std::size_t strpos,
            const std::size_t& length,
            const std::vector< std::string >& entropies
            )
    {
        double entropy = 0.0;

        if( strpos + length > entropies.size() )
            strpos = strpos - ( strpos + length - entropies.size() );

        for( std::size_t i = 0; i < length; i++ ) entropy += std::stod( entropies[ strpos + i ]);
        return entropy / (double)length;
    }
};

} // end of namespace algorithm
} // end of namespace ago
