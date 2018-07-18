#pragma once
#include <AGO/engine/components/named_component.hpp>
#include <CCD/para_thread_pool/para_thread_pool.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>
#include <mutex>

namespace ago {
namespace component {

class VcfAnnotationAnalysis : public engine::NamedComponent
{
    using Base = engine::NamedComponent;
    using VcfType = std::tuple< std::string, std::size_t, std::string, std::string, double, std::string >;
    //                              chr         pos         ref             alt     qul         id

    struct AnnotationVcf : public RawBedBase
    {
        char strand_;
        uint64_t end_;
        size_t is_filtered_;

        std::string chromosome_;
        std::vector< std::vector< std::string >> annotation_info_;

        std::string gene;
        std::string rmsk;
        std::string dbsnp;
        std::string comonsnp;
        std::string pub;
        std::string omim;
        std::string clinvar;
        std::string genedis;
        std::string cosmic;
        std::string struvart;
        std::string exac;
        std::string esp;
        std::string kgneome;

        VcfType vcf;

    	AnnotationVcf( void )
            : RawBedBase()
        {}
        
    	AnnotationVcf( const VcfType& in )
            : RawBedBase( std::get<0>( in ), std::get<1>( in ) -1, std::get<1>( in ) -1 + std::get<2>( in ).length(), "+" )
    		, strand_ ( '+' )
    		, end_( std::get<1>( in ) -1 + std::get<2>( in ).length() )
    		, chromosome_( std::get<0>( in ))
    		, annotation_info_(0)
    		, is_filtered_(0)
            , gene    (".")
            , rmsk    (".")
            , dbsnp   (".")
            , comonsnp("N")
            , pub     ("N")
            , omim    (".")
            , clinvar (".")
            , genedis (".")
            , cosmic  (".")
            , struvart(".")
            , exac    (".")
            , esp     ("N")
            , kgneome ("N")
            , vcf( in )
    	{}

        friend std::ostream& operator<< (std::ostream& out, const AnnotationVcf& vcf )
        {
            out << std::get<0>( vcf.vcf ) << "\t"
                << std::get<1>( vcf.vcf ) << "\t"
                << std::get<2>( vcf.vcf ) << "\t"
                << std::get<3>( vcf.vcf ) << "\t"
                << std::get<4>( vcf.vcf ) << "\t"
                << vcf.gene     << "\t"
                << vcf.rmsk     << "\t"
                << vcf.dbsnp    << "\t"
                << vcf.comonsnp << "\t"
                << vcf.pub      << "\t"
                << vcf.omim     << "\t"
                << vcf.clinvar  << "\t"
                << vcf.genedis  << "\t"
                << vcf.cosmic   << "\t"
                << vcf.struvart << "\t"
                << vcf.exac     << "\t"
                << vcf.esp      << "\t"
                << vcf.kgneome  << "\n";
            return out;
        }
    };

    struct HitHandler
    {
        template< class DB_BED, class ANN_BED >
        static void run( DB_BED&& db_bed, ANN_BED&& ann_bed, int& db_idx )
        {
            auto& type = std::get<4>( db_bed.data );
            auto& name = std::get<5>( db_bed.data );

            if( type == "Gencode" ) ann_bed.gene = name;
            if( type == "RMSK" || type == "NonRMSKofDuplications1000Base" ) ann_bed.rmsk = "Y";
            if( type == "Snp150All" ) ann_bed.dbsnp = name;
            if( type == "Snp150Common" ) ann_bed.comonsnp = "Y";
            if( type == "PubsMarkerSnp" ) ann_bed.pub = "Y";
            if( type == "Omim" ) ann_bed.omim = name;
            if( type == "Clinvar" ) ann_bed.clinvar = name;
            if( type == "GeneticDisorders" ) ann_bed.genedis = name;
            if( type == "COSMIC" ) ann_bed.cosmic = name;
            if( type == "StructuralVariation" ) ann_bed.struvart = name;
            if( type == "ExAC" ) ann_bed.exac = "Y";
            if( type == "ExonSequencingProject" ) ann_bed.esp = name;
            if( type == "1kGenome" ) ann_bed.kgneome = "Y";
        }
    };

    using BedFileReaderImpl =  FileReader_impl<
          Bed
        , std::tuple< std::string, uint32_t, uint32_t, char, std::string, std::string >
        , SOURCE_TYPE::IFSTREAM_TYPE
    >;

    using AnnotationTrait = Annotation<
          BedFileReaderImpl
        , AnnoIgnoreStrand::IGNORE
        , AnnoType::INTERSET
        , HitHandler
    >;

    using Annotations = AnnotationSet <
          AnnotationVcf
        , AnnotationTrait
    >;

    std::vector< std::string > annotation_files_;

    int task_number_;
    int thread_num_;

  protected:

    virtual void config_parameters( const bpt::ptree& p ) override
    {
        auto& db( this->mut_data_pool() );
        auto& pipeline_schema (db.pipeline_schema() );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "sample_files" ))
            db.push_path( "sample_files", child.second );

        for( auto& child : pipeline_schema.get_child( "input" ).get_child( "annotation_files" ))
            db.push_path( "annotation_files", child.second );

        task_number_ = p.get_optional< int >( "task_number" ).value_or( 50000 );
        thread_num_  = p.get_optional< int >( "thread_num" ).value_or( 16 );
    }

  public:

    using Base::Base;

    virtual void initialize() override
    {
        auto& db( this->mut_data_pool() );
        db.require_genome( db );

        if( db.exist_path_tag( "annotation_files" ))
        {
            for( auto& bed : db.get_path_list( "annotation_files" ))
            {
                annotation_files_.emplace_back( bed.string() );
            }
        }
    }

    virtual void start() override
    {
        auto& db( this->mut_data_pool() );
        auto& monitor = db.monitor();

        Annotations annotator( annotation_files_ );

        std::vector< std::string > vcf_paths( get_path_list_string( db.get_path_list( "sample_files" )));
        std::vector< std::vector< AnnotationVcf >> vcf_smp_pool( vcf_paths.size(), std::vector< AnnotationVcf >() );

        std::vector< std::pair< size_t, AnnotationVcf >> vcf_pool;
        std::map< std::thread::id, std::vector< std::pair< size_t, AnnotationVcf >>> vcf_tread_pool;

        monitor.set_monitor( "Component VcfAnnotationAnalysis", 5 );
        monitor.log( "Component VcfAnnotationAnalysis", "Reading Data" );

        ParaThreadPool smp_parallel_pool( vcf_paths.size() );
        ParaThreadPool vcf_parallel_pool( thread_num_ );
        std::mutex vcf_mutex;

        for( size_t smp = 0; smp < vcf_paths.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &vcf_paths, &vcf_smp_pool, this ] ()
            {
                std::string line;
                std::vector< std::string > splits;
                std::vector< std::string > chr_size;

                std::string sample_name = get_sample_name( vcf_paths[ smp ]);
                std::fstream filein( vcf_paths[ smp ], std::ios::in );

                while( std::getline( filein, line ))
                {
                    if( line.substr( 0, 1 ) == "#" ) continue;
                    boost::iter_split( splits, line, boost::algorithm::first_finder( "\t" ));
                    boost::iter_split( chr_size, splits[0], boost::algorithm::first_finder( "_" ));
                    if( chr_size.size() != 1 ) continue;

                    vcf_smp_pool[ smp ].emplace_back( AnnotationVcf( VcfType{
                            ( splits[0].length() > 3 ? splits[0] : "chr" + splits[0] )
                            , std::stol( splits[1] )
                            , splits[3]
                            , splits[4]
                            , std::stod( splits[5] )
                            , splits[2]
                            }));
                }

                filein.close();
            });
        }

        smp_parallel_pool.flush_pool();
        monitor.log( "Component VcfAnnotationAnalysis", "Annotating Data" );

        for( size_t smp = 0; smp < vcf_paths.size(); ++smp )
        {
            for( auto& vcf : vcf_smp_pool[ smp ] )
            {
                vcf_pool.emplace_back( std::make_pair( smp, vcf ));

                if( vcf_pool.size() >= task_number_ )
                {
                    vcf_parallel_pool.job_post([ vcf_pool, &vcf_tread_pool, &annotator, &vcf_mutex, this ] () mutable
                    {
                        {
                            std::lock_guard< std::mutex > smp_lock( vcf_mutex );
                            if( vcf_tread_pool.find( std::this_thread::get_id() ) == vcf_tread_pool.end() )
                                vcf_tread_pool[ std::this_thread::get_id() ] = std::vector< std::pair< size_t, AnnotationVcf >>(); 
                        }

                        for( auto& anno_vcf : vcf_pool )
                        {
                            annotator.AnnotateAll( anno_vcf.second );
                            vcf_tread_pool[ std::this_thread::get_id() ].emplace_back(
                                std::make_pair( anno_vcf.first, anno_vcf.second ));
                        }
                    });

                    vcf_pool.clear();
                }
            }
        }

        vcf_parallel_pool.job_post([ vcf_pool, &vcf_tread_pool, &annotator, &vcf_mutex, this ] () mutable
        {
            {
                std::lock_guard< std::mutex > vcf_lock( vcf_mutex );
                if( vcf_tread_pool.find( std::this_thread::get_id() ) == vcf_tread_pool.end() )
                    vcf_tread_pool[ std::this_thread::get_id() ] = std::vector< std::pair< size_t, AnnotationVcf >>(); 
            }

            for( auto& anno_vcf : vcf_pool )
            {
                annotator.AnnotateAll( anno_vcf.second );
                vcf_tread_pool[ std::this_thread::get_id() ].emplace_back(
                    std::make_pair( anno_vcf.first, anno_vcf.second ));
            }
        });

        vcf_pool.clear();
        vcf_smp_pool = std::vector< std::vector< AnnotationVcf >>( vcf_paths.size(), std::vector< AnnotationVcf >() );

        vcf_parallel_pool.flush_pool();
        monitor.log( "Component VcfAnnotationAnalysis", "Transforming Data" );

        for( auto& thread : vcf_tread_pool )
        {
            for( auto& smp_vcf : thread.second )
            {
                vcf_smp_pool[ smp_vcf.first ].emplace_back( smp_vcf.second );
            }
        }

        monitor.log( "Component VcfAnnotationAnalysis", "Outputing Data" );

        for( size_t smp = 0; smp < vcf_paths.size(); ++smp )
        {
            smp_parallel_pool.job_post([ smp, &vcf_paths, &vcf_smp_pool, this ] ()
            {
                std::ofstream output( get_sample_name( vcf_paths[ smp ] ) + ".text" );
                output << "Chr\tPos\tRef\tAlt\tQul\tGene\tRmsk\tdbSNP\tCommonSNP\tPublication\tOMIM\tClinvar\tGeneticDisorders\tCosmic\tStructuralVariation\tExAC\tExonSequencingProject\t1000Genome\n";
                for( auto& anno_vcf : vcf_smp_pool[ smp ]) output << anno_vcf;
                output.close();
            });
        }

        smp_parallel_pool.flush_pool();
        monitor.log( "Component VcfAnnotationAnalysis", "Complete" );
    }

    std::vector< std::string > get_path_list_string( const std::vector< boost::filesystem::path >& paths )
    {
        std::vector< std::string > res;
        for( auto& path : paths ) res.emplace_back( path.string() );
        return res;
    }

    std::string get_sample_name( const std::string& path )
    {
        std::vector< std::string > path_file;
        boost::iter_split( path_file, path, boost::algorithm::first_finder( "/" ));

        std::vector< std::string > sample;
        boost::iter_split( sample, path_file[ path_file.size()-1 ], boost::algorithm::first_finder( "." ));

        return sample[0];
    }
};

} // end of namespace component
} // end of namespace ago
