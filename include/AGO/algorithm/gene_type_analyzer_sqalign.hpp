#pragma once
#include <unistd.h>
#include <AGO/format/md_rawbed.hpp>

namespace ago {
namespace algorithm {

struct SeqType
{
    //                                  seed        arm         tail
    using SeedTailType = std::tuple< std::string, std::string, std::string >;
    using ReadType = std::tuple< std::size_t, std::size_t, double, char, std::size_t, std::map< std::size_t, char >, std::set< std::size_t >>;
    //                              start       length      ppm   isfilter  index       md_map    index       nt        tc_set      index
    
    char strand;

    std::string biotype;
    std::string fullseq;
    std::string chr;

    std::size_t refstart;
    std::size_t start;
    std::size_t end;

    std::map< std::string, std::string > *genome;
    std::map< SeedTailType, std::size_t > seed2read_idx;
    std::map< std::size_t, SeedTailType > read2seed_idx;

    std::vector< ReadType > reads_vec;

    SeqType()
    {}

    // void init( AnnotationRawBed<>& raw_bed, auto& genome_table, const std::string& biotype_ )
    void init( ago::format::MDRawBed& raw_bed, auto& genome_table, const std::string& biotype_ )
    {
        this->chr     = raw_bed.chromosome_;
        this->start   = raw_bed.start_;
        this->end     = raw_bed.end_;
        this->strand  = raw_bed.strand_;
        this->genome  = &genome_table;
        this->fullseq = "";
        this->biotype = biotype_;
        this->seed2read_idx.clear();
        this->read2seed_idx.clear();
        this->reads_vec.clear();
    }

    // std::string get_arm( AnnotationRawBed<>& raw_bed )
    std::string get_arm( ago::format::MDRawBed& raw_bed )
    {
        std::string arm;
        // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
        {
            std::size_t i = 0; // do first priority
            if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
            {
                if(( raw_bed.annotation_info_[i][0] == biotype ) ||
                   ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))) 
                {
                    for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                    {
                        arm = biotype.substr( 0, 5 ) == "miRNA" || biotype == "mirtron"
                            ? raw_bed.annotation_info_[i][ j+1 ].substr( raw_bed.annotation_info_[i][ j+1 ].length() -2, 2 )
                            : "." ;
                    }
                }
            }
        }
        return arm;
    }

    // SeedTailType make_seed_tail_index( AnnotationRawBed<>& rawbed )
    SeedTailType make_seed_tail_index( ago::format::MDRawBed& rawbed )
    {
        SeedTailType seedtemp = {
            GeneTypeAnalyzerCounting::seqT2U( rawbed.getReadSeq( *genome ).substr( 1, 7 )),
            get_arm( rawbed ),
            ( GeneTypeAnalyzerCounting::seqT2U( rawbed.getTail() ) != "" ? GeneTypeAnalyzerCounting::seqT2U( rawbed.getTail() ) : "." )
        };
        return seedtemp;
    }

    // void insert( AnnotationRawBed<>& rawbed, const double& ppm )
    void insert( ago::format::MDRawBed& rawbed )
    {
        if( rawbed.end_   > this->end   ) this->end   = rawbed.end_;
        if( rawbed.start_ < this->start ) this->start = rawbed.start_;

        SeedTailType seedtemp = make_seed_tail_index( rawbed );

        ReadType readtemp = {
            ( this->strand == '+' ? rawbed.start_ : rawbed.end_ ), 
            (int)rawbed.length_ - (int)rawbed.tail_length_,
            rawbed.ppm_,
            ( rawbed.is_filtered_ == 0 ? 'N' : 'Y' ),
            reads_vec.size(),
            rawbed.md_map,
            rawbed.tc_set
        };

        seed2read_idx[ seedtemp ] = reads_vec.size();
        read2seed_idx[ reads_vec.size() ] = seedtemp;

        reads_vec.emplace_back( readtemp );
    }

    void sorting_reads()
    {
        std::sort( reads_vec.begin(), reads_vec.end(), [ this ]( const ReadType& a, const ReadType& b )
        {
            if( std::get<0>(a) == std::get<0>(b) )
                if( std::get<1>(a) == std::get<1>(b) )
                    if( std::get<2>(a) == std::get<2>(b) )
                         return std::get<4>(a) < std::get<4>(b);
                    else return std::get<2>(a) > std::get<2>(b);
                else return std::get<1>(a) > std::get<1>(b);
            else switch( this->strand )
            {
                case '+' : return std::get<0>(a) < std::get<0>(b);
                case '-' : return std::get<0>(a) > std::get<0>(b);
            }
        });
    }

    void get_sequence()
    {
        std::string read_seq = ( this->strand == '+'
                ? genome->operator[]( this->chr ).substr( this->start -1, this->end - this->start )
                : genome->operator[]( this->chr ).substr( this->start -1, this->end - this->start )
                );
        std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper );

        if( this->strand == '-' )
        {
            std::transform( read_seq.begin(), read_seq.end(), read_seq.begin(), [ this ]( char c )
            {
                switch (c) {
                    case 'A': c = 'U'; break;
                    case 'U': c = 'A'; break;
                    case 'C': c = 'G'; break;
                    case 'G': c = 'C'; break;
                    case 'N': c = 'N'; break;
                }
                return c;
            });

            std::reverse( read_seq.begin(), read_seq.end() );
        }

        this->fullseq = read_seq;
    }

    void formation( std::size_t extend_refseq )
    {
        sorting_reads();
        std::map< std::size_t, SeedTailType > read2seed_idx_temp;

        for( std::size_t i = 0; i < reads_vec.size(); ++i )
        {
            seed2read_idx[ read2seed_idx[ std::get<4>( reads_vec[i] ) ]] = i;
            read2seed_idx_temp[ i ] = read2seed_idx[ std::get<4>( reads_vec[i] ) ];
            std::get<4>( reads_vec[i] ) = i;
            std::get<0>( reads_vec[i] ) = ( this->strand == '+'
                    ? std::get<0>( reads_vec[i] ) - this->start
                    : this->end - std::get<0>( reads_vec[i] )
                    ) + extend_refseq;
        }

        this->start = this->start - extend_refseq;
        this->end   = this->end   + extend_refseq;

        this->refstart = ( this->strand == '+' ? this->start : this->end );

        read2seed_idx = read2seed_idx_temp;
        get_sequence();
    }

    friend std::ostream& operator<< ( std::ostream& out, SeqType& seqs )
    {
        out << "\t" << seqs.fullseq << "\t" << seqs.chr << "\t" << seqs.refstart << "\t" << seqs.strand << "\n";
        for( auto& seq : seqs.reads_vec )
        {
            std::string md_tag = "";
            std::string tc_tag = "";

            if( std::get<5>( seq ).size() != 0 )
            {
                for( auto& md : std::get<5>( seq ))
                    md_tag += std::to_string( md.first ) + md.second;
            }

            if( std::get<6>( seq ).size() != 0 )
            {
                for( auto& tc : std::get<6>( seq ))
                    tc_tag += std::to_string( tc ) + 'C';
            }

            out << "\t" << std::get<0>( seq )   // start
                << "\t" << std::get<1>( seq )   // length
                << "\t" << std::get<2>( seq )   // ppm
                << "\t" << std::get<3>( seq )   // isfilter
                // << "\t" << std::get<0>( seqs.read2seed_idx[ std::get<4>( seq ) ])    // seed
                << "\t" << std::get<1>( seqs.read2seed_idx[ std::get<4>( seq ) ])       // arm
                << "\t" << std::get<2>( seqs.read2seed_idx[ std::get<4>( seq ) ])       // tail
                << "\t" << ( std::get<5>( seq ).size() == 0 ? "." : md_tag )            // MDtag
                << "\t" << ( std::get<6>( seq ).size() == 0 ? "." : tc_tag )            // TCtag
                << "\n";
        }

        return out;
    }

    std::size_t line_counts()
    {
        return this->reads_vec.size() + 1;
    }
};

class GeneTypeAnalyzerSqalign
{

  public:

    GeneTypeAnalyzerSqalign()
    {}

    static void output_sqalign(
            const std::string& output_name,
            std::vector< BedSampleType >& bed_samples,
            const std::string& biotype,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::size_t& extend_refseq,
            const std::size_t& max_anno_merge_size,
            const std::string& rnafold_path
            )
    {
        std::ofstream outidx, outtsv;
        std::map< std::string, SeqType > chr_mapping;
        std::map< std::string, std::vector< std::string >> rnafolds;

        std::size_t idx;
        double ppm = 0.0;

        if(( biotype == "miRNA_mirtron" || biotype == "miRNA" || biotype == "mirtron" ) && rnafold_path != "" )
        {
            get_rnafold( rnafold_path, rnafolds );
            if( !boost::filesystem::exists( output_name + "rnafold.tsv" ))
                boost::filesystem::create_symlink( rnafold_path, ( output_name + "rnafold.tsv" ).c_str() );
        }

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            idx = 0;
            chr_mapping = get_chrmap_table( bed_samples[ smp ].second, biotype, genome_table, filter_ppm, extend_refseq, max_anno_merge_size );

            if(( biotype == "miRNA_mirtron" || biotype == "miRNA" || biotype == "mirtron" ) && rnafold_path != "" )
                apply_rnafold( chr_mapping, rnafolds, extend_refseq );

            outidx.open( output_name + bed_samples[ smp ].first + ".idx" );
            outtsv.open( output_name + bed_samples[ smp ].first + ".tsv" );

            for( auto& mir : chr_mapping )
            {
                ppm = get_region_ppm( mir.second );
                outidx << mir.first << "\t" << idx << "\t" << ppm << "\n";;
                outtsv << mir.first << mir.second;
                idx += mir.second.line_counts();
            }

            outidx.close();
            outtsv.close();
        }
    }

    static double get_region_ppm( SeqType& seqs )
    {
        double ppm = 0.0;
        for( auto& seq : seqs.reads_vec )
            ppm += std::get<2>( seq );
        return ppm;
    }

    static void get_rnafold( const std::string& rnafold_path, std::map< std::string, std::vector< std::string >>& rnafolds )
    {
        std::fstream file( rnafold_path, std::ios::in );
        std::vector< std::string > split;
        std::string line;

        while( std::getline( file, line ))
        {
            boost::iter_split( split, line, boost::algorithm::first_finder( "\t" ));
            rnafolds[ split[0] ] = split;
        }
    }

    static void apply_rnafold( std::map< std::string, SeqType >& chr_mapping, std::map< std::string, std::vector< std::string >>& rnafolds, const std::size_t& extend_refseq )
    {
        double shift;
        std::size_t start;

        for( auto& mir : chr_mapping )
        {
            if( rnafolds.find( mir.first ) == rnafolds.end() ) continue;

            start = std::stoll( mir.second.strand == '+' ? rnafolds[ mir.first ][1] : rnafolds[ mir.first ][2] ) + 1;
            shift = start > mir.second.refstart ? ( start - mir.second.refstart ) : ( mir.second.refstart - start );

            for( auto& seq : mir.second.reads_vec )
                std::get<0>( seq ) 
                    = mir.second.strand == '+'
                    ? ( start > mir.second.refstart ? ( std::get<0>( seq ) - shift ) : ( std::get<0>( seq ) + shift ))
                    : ( start > mir.second.refstart ? ( std::get<0>( seq ) + shift ) : ( std::get<0>( seq ) - shift ))
                    ;

            mir.second.refstart = start;
            mir.second.fullseq = rnafolds[ mir.first ][4];
        }
    }

    static std::map< std::string, SeqType > get_chrmap_table(
            // std::vector< AnnotationRawBed<> >& smp_anno,
            std::vector< ago::format::MDRawBed >& smp_anno,
            const std::string& biotype,
            auto& genome_table,
            const std::size_t& filter_ppm,
            const std::size_t& extend_refseq,
            const std::size_t& max_anno_merge_size
            )
    {
        std::map< std::string, std::vector< ago::format::MDRawBed >> anno_beds;
        std::map< std::string, SeqType > chr_mappings;
        std::string anno;

        for( auto& raw_bed : smp_anno )
        {
            if( raw_bed.ppm_ < filter_ppm ) continue;
            if( !raw_bed.annotation_info_.empty() && !raw_bed.annotation_info_[0].empty() && raw_bed.start_ > extend_refseq )
            {
                if( biotype != "miRNA" && biotype != "mirtron" && biotype != "miRNA_mirtron" && biotype == raw_bed.annotation_info_[0][0] )
                {
                    anno = raw_bed.annotation_info_[0][1];
                    anno_beds[ anno ].emplace_back( raw_bed );
                    continue;
                }

                if(( raw_bed.annotation_info_[0][0] == biotype ) ||
                   ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[0][0] == "miRNA" || raw_bed.annotation_info_[0][0] == "mirtron" ))) 
                {
                        anno = raw_bed.annotation_info_[0][1].substr( 0, raw_bed.annotation_info_[0][1].length() -3 );
                        if( chr_mappings.find( anno ) == chr_mappings.end() )
                        {
                            chr_mappings[ anno ] = SeqType();
                            chr_mappings[ anno ].init( raw_bed, genome_table, raw_bed.annotation_info_[0][0] );
                        }
                        chr_mappings[ anno ].insert( raw_bed );
                }
            }
        }

        for( auto& annobed : anno_beds ) annotation_reformating( annobed.second, max_anno_merge_size );
        for( auto& annobed : anno_beds ) for( auto& bed : annobed.second )
        {
            anno = bed.annotation_info_[0][1];
            if( chr_mappings.find( anno ) == chr_mappings.end() )
            {
                chr_mappings[ anno ] = SeqType();
                chr_mappings[ anno ].init( bed, genome_table, bed.annotation_info_[0][0] );
            }
            chr_mappings[ anno ].insert( bed );
        }

        for( auto& mir : chr_mappings ) mir.second.formation( extend_refseq );
        return chr_mappings;
    }

    static void annotation_reformating( std::vector< ago::format::MDRawBed >& beds, const std::size_t& max_anno_merge_size )
    {
        std::string region = "";
        std::set< std::tuple< std::string, std::size_t, std::size_t, std::string >> gene_region_set;
        std::map< std::string, std::size_t > gene_counts;

        //          regionID                    chr         start       end         anno
        std::map< std::string, std::tuple< std::string, std::size_t, std::size_t, std::string >> bed_list;
        std::map< std::string, std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>> bed_map;
        //          chr                     start                 end                       regionID
    
        for( auto& bed : beds )
        {
            region = bed.chromosome_ + ":" + std::to_string( bed.start_ ) + "-" + std::to_string( bed.end_ );

            if( bed_list.find( region ) == bed_list.end() )
            {
                bed_list[ region ] = std::make_tuple( bed.chromosome_, bed.start_, bed.end_, bed.annotation_info_[0][1] );
                bed_map[ bed.chromosome_ ][ bed.start_ ][ bed.end_ ].emplace_back( region );
            }
        }

        for( auto& bed : bed_map  ) formating_start_end( bed.second, bed_list, max_anno_merge_size );
        for( auto& bed : bed_list ) gene_region_set.emplace( bed.second );

        for( auto& bed : gene_region_set )
        {
            if( gene_counts.find( std::get<3>( bed )) == gene_counts.end() )
                gene_counts[ std::get<3>( bed ) ] = 0.0;
            gene_counts[ std::get<3>( bed ) ]++;
        }

        for( auto& bed : beds )
        {
            region = bed.chromosome_ + ":" + std::to_string( bed.start_ ) + "-" + std::to_string( bed.end_ );

            if( gene_counts[ bed.annotation_info_[0][1] ] != 1 )
            {
                bed.annotation_info_[0][1] =
                    bed.annotation_info_[0][1] + "|" +
                    std::get<0>( bed_list[ region ]) + ":" +
                    std::to_string( std::get<1>( bed_list[ region ])) + "-" +
                    std::to_string( std::get<2>( bed_list[ region ])) ;
            }
        }
    }

    static void formating_start_end(
            std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>& beds,
            std::map< std::string, std::tuple< std::string, std::size_t, std::size_t, std::string >>& bed_list,
            const std::size_t& max_anno_merge_size
            )
    {
        std::size_t start = beds.begin()->first;
        std::size_t end = beds.begin()->second.begin()->first;

        for( auto sit = beds.begin(); sit != beds.end(); ++sit )
            for( auto eit = sit->second.begin(); eit != sit->second.end(); ++eit )
            {
                if( eit->first - start > max_anno_merge_size && end < sit->first )
                {
                    rerange_start_end(
                            beds,
                            bed_list,
                            max_anno_merge_size,
                            start,
                            end
                            );
                    
                    start = sit->first;
                }

                end = eit->first;
            }

        rerange_start_end(
                beds,
                bed_list,
                max_anno_merge_size,
                start,
                end
                );
    }

    static void rerange_start_end(
            std::map< std::size_t, std::map< std::size_t, std::vector< std::string >>>& beds,
            std::map< std::string, std::tuple< std::string, std::size_t, std::size_t, std::string >>& bed_list,
            const std::size_t& max_anno_merge_size,
            const std::size_t& start,
            const std::size_t& end
            )
    {
        for( auto sit = beds.find( start ); sit != beds.end(); ++sit )
            for( auto eit = sit->second.begin(); eit != sit->second.end(); ++eit )
            {
                if( end - start > max_anno_merge_size && end < sit->first ) return;
                for( auto& region_id : eit->second )
                {
                    std::get<1>( bed_list[ region_id ] ) = start;
                    std::get<2>( bed_list[ region_id ] ) = end;
                }
            }
    }

    static void output_sqalign_visualization( const std::string& output_name, const std::size_t& extend_refseq )
    {
        std::ofstream output( output_name + "index.php" );

        output << "<!DOCTYPE html>" << "\n";
        output << "<html>" << "\n";
        output << "    <meta charset='utf-8'>" << "\n";
        output << "    <script src='https://d3js.org/d3.v5.min.js'></script>" << "\n";
        output << "    <script src='https://code.jquery.com/jquery-3.3.1.min.js' ></script>" << "\n";
        output << "    <script src='https://code.jquery.com/ui/1.12.1/jquery-ui.js'></script>" << "\n";
        output << "    <link rel='stylesheet' href='https://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css'>" << "\n";
        output << "    <link href='https://fonts.googleapis.com/css?family=Source+Code+Pro' rel='stylesheet'>" << "\n";
        output << "    <style>.ui-tooltip-content{font-size:15px;font-family:Calibri;}</style>" << "\n";
        output << "    <body>" << "\n";
        output << "    <?" << "\n";
        output << "        Shell_Exec( 'rm /tmp/*' );" << "\n";
        output << "" << "\n";
        output << "        $TSV_File = $_GET['TSV_File'];" << "\n";
        output << "        $Annotation_Select = $_GET['Annotation_Select'];" << "\n";
        output << "        $Seed_Select = $_GET['Seed_Select'];" << "\n";
        output << "" << "\n";
        output << "        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        if( $_GET['Annotation_Select'] == '' ) $Annotation_Select = $_POST['Annotation_Select'];" << "\n";
        output << "        if( $_GET['Seed_Select'] == '' ) $Seed_Select = $_POST['Seed_Select'];" << "\n";
        output << "" << "\n";
        output << "        $SVG = $_POST['SVG'];" << "\n";
        output << "        $Length = $_POST['Length'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $MaxPPM = $_POST['MaxPPM'];" << "\n";
        output << "        $SumPPM = $_POST['SumPPM'];" << "\n";
        output << "        $Color_Hight = $_POST['Color_Hight'];" << "\n";
        output << "        $Color_Low = $_POST['Color_Low'];" << "\n";
        output << "        $Segment_Select = $_POST['Segment_Select'];" << "\n";
        output << "" << "\n";
        output << "        if( $Color_Hight == '' ) $Color_Hight = 'Black';" << "\n";
        output << "        if( $Color_Low == '' ) $Color_Low = 'WhiteSmoke';" << "\n";
        output << "" << "\n";
        output << "        $isSVG = $SVG == 'Yes' ? true : false ;" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "        $isRNAfold = false;" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            if( $TSV_List[$i] == 'rnafold.tsv' )" << "\n";
        output << "            {" << "\n";
        output << "                $isRNAfold = true;" << "\n";
        output << "                continue;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='SVG' value='$SVG' />" << "\n";
        output << "            <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='MaxPPM' value='$MaxPPM' />" << "\n";
        output << "            <input type='hidden' name='SumPPM' value='$SumPPM' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Seed_Select' value='$Seed_Select' />" << "\n";
        output << "            <input type='hidden' name='Segment_Select' value='$Segment_Select' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( substr( $TSV_File, 0, strlen( $TSV_File ) -4 ).'.idx' );" << "\n";
        output << "        $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
        output << "" << "\n";
        output << "        $Total_PPM = 0;" << "\n";
        output << "        $PPM_Array = Array();" << "\n";
        output << "        $Def_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $Anno_Array = Array();" << "\n";
        output << "        $Data_Array = Array();" << "\n";
        output << "" << "\n";
        output << "        $Ref_Chr    = '';" << "\n";
        output << "        $Ref_Start  = '';" << "\n";
        output << "        $Ref_Strand = '';" << "\n";
        output << "        $Full_Sequc = '';" << "\n";
        output << "        $Data_Check = false;" << "\n";
        output << "" << "\n";
        output << "        For( $j = 0; $j < Count( $inFile_Lines )-1; ++$j )" << "\n";
        output << "        {" << "\n";
        output << "            $anno = Explode( \"\\t\", $inFile_Lines[$j] );" << "\n";
        output << "            $Anno_Array[ $anno[0] ] = $anno[1];" << "\n";
        output << "" << "\n";
        output << "            $Anno_Split = Explode( '|', $anno[0] );" << "\n";
        output << "            if( Count( $Anno_Split ) != 1 )" << "\n";
        output << "            {" << "\n";
        output << "                if( Array_Key_Exists( $Anno_Split[0], $Def_Array ))" << "\n";
        output << "                {" << "\n";
        output << "                    if( $anno[2] > $PPM_Array[ $Def_Array[ $Anno_Split[0] ]] )" << "\n";
        output << "                        $Def_Array[ $Anno_Split[0] ] = $anno[0];" << "\n";
        output << "                }" << "\n";
        output << "                else $Def_Array[ $Anno_Split[0] ] = $anno[0];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $PPM_Array [ $anno[0] ] = $anno[2];" << "\n";
        output << "            if( $Total_PPM < $anno[2] ) $Total_PPM = $anno[2];" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Ksort( $Anno_Array );" << "\n";
        output << "        Ksort( $PPM_Array );" << "\n";
        output << "" << "\n";
        output << "        if( Array_Key_Exists( $Annotation_Select, $Def_Array ))" << "\n";
        output << "            $Annotation_Select = $Def_Array[ $Annotation_Select ];" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $idx )" << "\n";
        output << "        {" << "\n";
        output << "            if( $Annotation_Select == $anno && $Annotation_Select != '' )" << "\n";
        output << "            {" << "\n";
        output << "                $TSV = new SplFileObject( $TSV_File );" << "\n";
        output << "                $TSV->seek( $idx );" << "\n";
        output << "" << "\n";
        output << "                $annoSeq = Explode( \"\\t\", trim( $TSV->current(), \"\\n\" ));" << "\n";
        output << "" << "\n";
        output << "                For( $i = 0; $i < Strlen( $annoSeq[1] ); $i++ )" << "\n";
        output << "                    $Full_Sequc = $Full_Sequc.( $annoSeq[1][$i] == 'T' ? 'U' : $annoSeq[1][$i] );" << "\n";
        output << "" << "\n";
        output << "                $Ref_Chr    = $annoSeq[2];" << "\n";
        output << "                $Ref_Start  = $annoSeq[3];" << "\n";
        output << "                $Ref_Strand = $annoSeq[4];" << "\n";
        output << "" << "\n";
        output << "                while( true )" << "\n";
        output << "                {" << "\n";
        output << "                    $TSV->next();" << "\n";
        output << "                    $tsvData = Explode( \"\\t\", trim( $TSV->current(), \"\\n\" ));" << "\n";
        output << "" << "\n";
        output << "                    if( $tsvData[0] != '' ) break;" << "\n";
        output << "                    $T2U = '';" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 0; $i < Strlen( $tsvData[6] ); $i++ )" << "\n";
        output << "                        $T2U = $T2U.( $tsvData[6][$i] == 'T' ? 'U' : $tsvData[6][$i] );" << "\n";
        output << "" << "\n";
        output << "                    $tsvData[6] = $T2U;" << "\n";
        output << "                    $T2U = '';" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 0; $i < Strlen( $tsvData[7] ); $i++ )" << "\n";
        output << "                        $T2U = $T2U.( $tsvData[7][$i] == 'T' ? 'U' : $tsvData[7][$i] );" << "\n";
        output << "" << "\n";
        output << "                    $tsvData[7] = $T2U;" << "\n";
        output << "                    $T2U = '';" << "\n";
        output << "" << "\n";
        output << "                    For( $i = 0; $i < Strlen( $tsvData[8] ); $i++ )" << "\n";
        output << "                        $T2U = $T2U.( $tsvData[8][$i] == 'T' ? 'U' : $tsvData[8][$i] );" << "\n";
        output << "" << "\n";
        output << "                    $tsvData[8] = $T2U;" << "\n";
        output << "                    Array_Push( $Data_Array, $tsvData );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $TSV = null;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Annotation_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $idx )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$anno.' ';" << "\n";
        output << "            if( $Annotation_Select == $anno )  echo 'selected ';" << "\n";
        output << "            echo '>'.$anno.' ('.Number_Format( (float)$PPM_Array[$anno], 2, '.', '' ).'ppm)</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='SVG' value='$SVG' />" << "\n";
        output << "            <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='MaxPPM' value='$MaxPPM' />" << "\n";
        output << "            <input type='hidden' name='SumPPM' value='$SumPPM' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "        $Total_PPM = 100 / $Total_PPM;" << "\n";
        output << "        echo \"<script>var select_color_map = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);\";" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $idx )" << "\n";
        output << "        {" << "\n";
        output << "           echo \"$( 'option[value=\\\"$anno\\\"]' ).css(" << "\n";
        output << "           {" << "\n";
        output << "               'background-color': select_color_map( '\".($PPM_Array[$anno] * $Total_PPM ).\"' )," << "\n";
        output << "               'color': '\".(( $PPM_Array[$anno] * $Total_PPM ) > 25 ? 'white' : 'black'  ).\"'" << "\n";
        output << "           });\";" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo '</script>';" << "\n";
        output << "" << "\n";
        output << "#<!--=================== RNA-Fold Datas ======================-->" << "\n";
        output << "" << "\n";
        output << "        $RNAfold = Array();" << "\n";
        output << "" << "\n";
        output << "        if( $isRNAfold && $Annotation_Select != '' )" << "\n";
        output << "        {" << "\n";
        output << "            $inFile = new SplFileObject( 'rnafold.tsv' );" << "\n";
        output << "" << "\n";
        output << "            while( !$inFile->eof() )" << "\n";
        output << "            {" << "\n";
        output << "                $inFile_Lines = $inFile->fgets();" << "\n";
        output << "                $inFile_Line = Explode( \"\\t\", Rtrim( $inFile_Lines ));" << "\n";
        output << "                if( $inFile_Line[0] != $Annotation_Select ) continue;" << "\n";
        output << "" << "\n";
        output << "                $Entropy = '';" << "\n";
        output << "                $Entropy_Max = 0;" << "\n";
        output << "" << "\n";
        output << "                $Sequ = '';" << "\n";
        output << "                $Fold = $inFile_Line[5];" << "\n";
        output << "" << "\n";
        output << "                For( $i = 0; $i < Strlen( $inFile_Line[4] ); $i++ )" << "\n";
        output << "                    $Sequ = $Sequ.( $inFile_Line[4][$i] == 'T' ? 'U' : $inFile_Line[4][$i] );" << "\n";
        output << "" << "\n";
        output << "                $Fold1 = '';" << "\n";
        output << "                $Fold2 = '';" << "\n";
        output << "                $Fold3 = ' ';" << "\n";
        output << "                $Fold4 = '';" << "\n";
        output << "                $Fold5 = '';" << "\n";
        output << "" << "\n";
        output << "                $Ltmp = '';" << "\n";
        output << "                $Mtmp = '';" << "\n";
        output << "                $Rtmp = '';" << "\n";
        output << "" << "\n";
        output << "                $isLloop = false;" << "\n";
        output << "                $isRloop = false;" << "\n";
        output << "" << "\n";
        output << "                $L = 0;" << "\n";
        output << "                $R = Strlen( $Fold )-1;" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i <= 20; ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Fold[$i] == ')' )  $isLloop = true;" << "\n";
        output << "                    if( $Fold[$i] == '(' && $isLloop )" << "\n";
        output << "                    {" << "\n";
        output << "                        $L = $i;" << "\n";
        output << "                        $isLloop = false;" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                for( $i = Strlen( $Fold )-1; $i >= Strlen( $Fold ) - 20; --$i )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Fold[$i] == '(' )  $isRloop = true;" << "\n";
        output << "                    if( $Fold[$i] == ')' && $isRloop )" << "\n";
        output << "                    {" << "\n";
        output << "                        $R = $i;" << "\n";
        output << "                        $isRloop = false;" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i < $L; ++$i )" << "\n";
        output << "                    $Ltmp = $Ltmp.$Sequ[$i];" << "\n";
        output << "" << "\n";
        output << "                for( $i = Strlen( $Fold )-1; $i > $R; --$i )" << "\n";
        output << "                    $Rtmp = $Rtmp.$Sequ[$i];" << "\n";
        output << "" << "\n";
        output << "                $xEntropy_Array = Array();" << "\n";
        output << "                $yEntropy_Array = Array();" << "\n";
        output << "                $vEntropy_Array = Array_Slice( $inFile_Line, 6 );" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i < Count( $vEntropy_Array ); ++$i )" << "\n";
        output << "                    if( $vEntropy_Array[$i] > $Entropy_Max ) $Entropy_Max = $vEntropy_Array[$i];" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i < Count( $vEntropy_Array ); ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $isSVG )" << "\n";
        output << "                    {" << "\n";
        output << "                        $xEntropy_Array[$i] = 9 + ( $i * 9 );" << "\n";
        output << "                        $yEntropy_Array[$i] = 46 - ( $vEntropy_Array[$i] * 46 / $Entropy_Max );" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $xEntropy_Array[$i] = 9 + ( $i * 18 );" << "\n";
        output << "                        $yEntropy_Array[$i] = 100 - ( $vEntropy_Array[$i] * 100 / $Entropy_Max );" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $Entropy = $Entropy.$xEntropy_Array[$i].','.$yEntropy_Array[$i].' ';" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i < Strlen( $Fold ); ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Fold[$L] == '(' && $Fold[$R] == ')' )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Spce = '';" << "\n";
        output << "                        $Lspc = '';" << "\n";
        output << "                        $Rspc = '';" << "\n";
        output << "" << "\n";
        output << "                        $Lcnt = Strlen( $Ltmp );" << "\n";
        output << "                        $Rcnt = Strlen( $Rtmp );" << "\n";
        output << "                        $Lnth = $Lcnt < $Rcnt ? $Rcnt : $Lcnt;" << "\n";
        output << "" << "\n";
        output << "                        for( $s = 0; $s < $Lnth; ++$s ) $Spce = $Spce.' ';" << "\n";
        output << "                        for( $s = 0; $s < $Rcnt - $Lcnt; ++$s ) $Lspc = $Lspc.' ';" << "\n";
        output << "                        for( $s = 0; $s < $Lcnt - $Rcnt; ++$s ) $Rspc = $Rspc.' ';" << "\n";
        output << "" << "\n";
        output << "                        $Fold1 = $Fold1.$Lspc.$Ltmp.' ';" << "\n";
        output << "                        $Fold2 = $Fold2.$Spce.$Sequ[$L];" << "\n";
        output << "" << "\n";
        output << "                        if(( $Sequ[$L] == 'A' && $Sequ[$R] == 'U' ) ||" << "\n";
        output << "                           ( $Sequ[$R] == 'A' && $Sequ[$L] == 'U' ) ||" << "\n";
        output << "                           ( $Sequ[$L] == 'C' && $Sequ[$R] == 'G' ) ||" << "\n";
        output << "                           ( $Sequ[$R] == 'C' && $Sequ[$L] == 'G' ) )" << "\n";
        output << "                           $Fold3 = $Fold3.$Spce.'|';" << "\n";
        output << "                        else" << "\n";
        output << "                           $Fold3 = $Fold3.$Spce.':';" << "\n";
        output << "" << "\n";
        output << "                        $Fold4 = $Fold4.$Spce.$Sequ[$R];" << "\n";
        output << "                        $Fold5 = $Fold5.$Rspc.$Rtmp.' ';" << "\n";
        output << "" << "\n";
        output << "                        $Ltmp = '';" << "\n";
        output << "                        $Rtmp = '';" << "\n";
        output << "" << "\n";
        output << "                        $L++;" << "\n";
        output << "                        $R--;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( $L == $R )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Mtmp = $Sequ[$L];" << "\n";
        output << "                        break;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( $Fold[$L] != '(' || $Fold[$R] != ')' )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( $Fold[$L] != '(' )" << "\n";
        output << "                        {" << "\n";
        output << "                            $Ltmp = $Ltmp.$Sequ[$L];" << "\n";
        output << "                            $L++;" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        if( $Fold[$R] != ')' )" << "\n";
        output << "                        {" << "\n";
        output << "                            $Rtmp = $Rtmp.$Sequ[$R];" << "\n";
        output << "                            $R--;" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( $L == ( $R + 1 )) break;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( Strlen( $Ltmp ) != Strlen( $Rtmp ))" << "\n";
        output << "                {" << "\n";
        output << "                    $Lstr = &$Ltmp;" << "\n";
        output << "                    $Sstr = &$Rtmp;" << "\n";
        output << "" << "\n";
        output << "                    if( Strlen( $Ltmp ) < Strlen( $Rtmp ))" << "\n";
        output << "                    {" << "\n";
        output << "                        $Lstr = &$Rtmp;" << "\n";
        output << "                        $Sstr = &$Ltmp;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $Dist = Strlen( $Lstr ) - Strlen( $Sstr );" << "\n";
        output << "" << "\n";
        output << "                    if( $Dist % 2 == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));" << "\n";
        output << "                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));" << "\n";
        output << "                    }" << "\n";
        output << "                    else if( $Mtmp != '' )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Sstr = $Sstr.$Mtmp;" << "\n";
        output << "                        $Dist = Strlen( $Lstr ) - Strlen( $Sstr );" << "\n";
        output << "                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));" << "\n";
        output << "                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));" << "\n";
        output << "                        $Mtmp = '';" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $Dist = $Dist -1;" << "\n";
        output << "                        $Sstr = $Sstr.Strrev( Substr( $Lstr, ( Strlen( $Lstr ) - ( $Dist / 2 ))));" << "\n";
        output << "                        $Mtmp = Substr( $Lstr, Strlen( $Sstr ));" << "\n";
        output << "                        $Lstr = Substr( $Lstr, 0, Strlen( $Sstr ));" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                for( $i = 0; $i < Strlen( $Ltmp ); ++$i )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $i < Strlen( $Ltmp ) -1 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $Fold1 = $Fold1.$Ltmp[$i];" << "\n";
        output << "                        $Fold2 = $Fold2.' ';" << "\n";
        output << "                        $Fold3 = $Fold3.' ';" << "\n";
        output << "                        $Fold4 = $Fold4.' ';" << "\n";
        output << "                        $Fold5 = $Fold5.$Rtmp[$i];" << "\n";
        output << "                    }" << "\n";
        output << "                    else" << "\n";
        output << "                    {" << "\n";
        output << "                        $Fold1 = $Fold1.' ';" << "\n";
        output << "                        $Fold2 = $Fold2.$Ltmp[$i];" << "\n";
        output << "                        $Fold3 = $Fold3.( $Mtmp == '' ? ' )' : ( ' '.$Mtmp ));" << "\n";
        output << "                        $Fold4 = $Fold4.$Rtmp[$i];" << "\n";
        output << "                        $Fold5 = $Fold5.' ';" << "\n";
        output << "                    }" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $RNAfold[ 'MFE' ][0] = $inFile_Line[3];" << "\n";
        output << "                $RNAfold[ 'Fold0' ][0] = $Fold;" << "\n";
        output << "                $RNAfold[ 'Fold1' ][0] = $Fold1;" << "\n";
        output << "                $RNAfold[ 'Fold2' ][0] = $Fold2;" << "\n";
        output << "                $RNAfold[ 'Fold3' ][0] = $Fold3;" << "\n";
        output << "                $RNAfold[ 'Fold4' ][0] = $Fold4;" << "\n";
        output << "                $RNAfold[ 'Fold5' ][0] = $Fold5;" << "\n";
        output << "                $RNAfold[ 'MaxEntropy' ][0] = $Entropy_Max;" << "\n";
        output << "                $RNAfold[ 'Entropy' ][0] = $Entropy;" << "\n";
        output << "                $RNAfold[ 'vEntropy' ] = $vEntropy_Array;" << "\n";
        output << "                $RNAfold[ 'xEntropy' ] = $xEntropy_Array;" << "\n";
        output << "                $RNAfold[ 'yEntropy' ] = $yEntropy_Array;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Select Segment ======================-->" << "\n";
        output << "" << "\n";
        output << "        $Seg_Start = 0;" << "\n";
        output << "        $Seg_End   = 0;" << "\n";
        output << "" << "\n";
        output << "        if( $Full_Sequc != '' && $Data_Array[0][5] == '.' )" << "\n";
        output << "        {" << "\n";
        output << "            $Segm_Uniqs = Array();" << "\n";
        output << "            $Segm_Array = Array();" << "\n";
        output << "            $Segm_Range = Array();" << "\n";
        output << "" << "\n";
        output << "            $Segm_Range[0] = 0; // start" << "\n";
        output << "            $Segm_Range[1] = 0; // end" << "\n";
        output << "            $Segm_Range[2] = 0; // ppm" << "\n";
        output << "" << "\n";
        output << "            $Total_PPM = 0;" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Data_Array ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                $End = $Data_Array[$i][1] + $Data_Array[$i][2]" << "\n";
        output << "                    +( $Data_Array[$i][6] != '.' ? Strlen( $Data_Array[$i][6] ) : 0 );" << "\n";
        output << "" << "\n";
        output << "                if( $Data_Array[$i][1] <= $Segm_Range[1] )" << "\n";
        output << "                {" << "\n";
        output << "                    if( $End > $Segm_Range[1] ) $Segm_Range[1] = $End;" << "\n";
        output << "                }" << "\n";
        output << "                else" << "\n";
        output << "                {" << "\n";
        output << "                    if( $Segm_Range[0] != $Segm_Range[1] )" << "\n";
        output << "                        Array_push( $Segm_Uniqs, $Segm_Range );" << "\n";
        output << "" << "\n";
        output << "                    for( $j = Count( $Segm_Array ); $j < $i; ++$j )" << "\n";
        output << "                        $Segm_Array[$j] = $Segm_Range;" << "\n";
        output << "" << "\n";
        output << "                    $Segm_Range[0] = $Data_Array[$i][1];" << "\n";
        output << "                    $Segm_Range[1] = $End;" << "\n";
        output << "                    $Segm_Range[2] = 0;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $Segm_Range[2] += $Data_Array[$i][3];" << "\n";
        output << "                if( $Total_PPM < $Data_Array[$i][3] )" << "\n";
        output << "                    $Total_PPM = $Data_Array[$i][3];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            $Total_PPM = 100 / $Total_PPM;" << "\n";
        output << "" << "\n";
        output << "            Array_push( $Segm_Uniqs, $Segm_Range );" << "\n";
        output << "" << "\n";
        output << "            if( $Segment_Select == '' )" << "\n";
        output << "            {" << "\n";
        output << "                $segMax = [ 0, 0, 0 ];" << "\n";
        output << "                for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )" << "\n";
        output << "                    if( $segMax[2] < $Segm_Uniqs[$i][2] ) $segMax = $Segm_Uniqs[$i];" << "\n";
        output << "                $Segment_Select = $segMax[0].'-'.$segMax[1];" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            for( $j = Count( $Segm_Array ); $j < Count( $Data_Array ); ++$j )" << "\n";
        output << "                $Segm_Array[$j] = $Segm_Range;" << "\n";
        output << "" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<select name=Segment_Select onchange=this.form.submit();>';" << "\n";
        output << "            echo \"<option value='' \"; if($Segment_Select=='') echo 'selected'; echo '>Select Segments</option>';" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value=\"'.$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].'\" ';" << "\n";
        output << "" << "\n";
        output << "                if( $Segment_Select == $Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1] )" << "\n";
        output << "                {" << "\n";
        output << "                    echo 'selected ';" << "\n";
        output << "                    $Seg_Start = $Segm_Uniqs[$i][0];" << "\n";
        output << "                    $Seg_End   = $Segm_Uniqs[$i][1];" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                echo '>'" << "\n";
        output << "                    .$Ref_Chr.':'" << "\n";
        output << "                    .( $Ref_Strand == '+'" << "\n";
        output << "                        ? (( (int)$Ref_Start + (int)$Segm_Uniqs[$i][0] ).'-'.( (int)$Ref_Start + (int)$Segm_Uniqs[$i][1] ))" << "\n";
        output << "                        : (( (int)$Ref_Start - (int)$Segm_Uniqs[$i][1] ).'-'.( (int)$Ref_Start - (int)$Segm_Uniqs[$i][0] ))" << "\n";
        output << "                    ).' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).'ppm)</option>';" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>" << "\n";
        output << "                <input type='hidden' name='SVG' value='$SVG' />" << "\n";
        output << "                <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "                <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "                <input type='hidden' name='MaxPPM' value='$MaxPPM' />" << "\n";
        output << "                <input type='hidden' name='SumPPM' value='$SumPPM' />" << "\n";
        output << "                <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "                <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "                <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "                <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "                <input type='hidden' name='Seed_Select' value='$Seed_Select' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>var select_color_map_seg = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);\";" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )" << "\n";
        output << "                echo \"$( 'option[value=\\\"\".$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).\"ppm)\\\"]' ).css(" << "\n";
        output << "                {" << "\n";
        output << "                    'background-color': select_color_map_seg( '\".( $Segm_Uniqs[$i][2] * $Total_PPM ).\"' )," << "\n";
        output << "                    'color': '\".(( $Segm_Uniqs[$i][2] * $Total_PPM ) > 25 ? 'white' : 'black'  ).\"'" << "\n";
        output << "                });\";" << "\n";
        output << "" << "\n";
        output << "            echo '</script>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Specific Seed ======================-->" << "\n";
        output << "" << "\n";
        output << "        $Seed_List = Array();" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Data_Array ); ++$i )" << "\n";
        output << "            Array_Push( $Seed_List, ( $Data_Array[$i][1] +2 ).':'.Substr( $Full_Sequc, $Data_Array[$i][1] +1, 7 ));" << "\n";
        output << "" << "\n";
        output << "        $Seed_List = Array_Unique( $Seed_List );" << "\n";
        output << "        $Seed_List = Array_Values( $Seed_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Seed_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Seed_Select=='') echo 'selected'; echo '>All Seeds</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < Count( $Seed_List ); ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$Seed_List[$i].' ';" << "\n";
        output << "            if( $Seed_Select == $Seed_List[$i] ) echo 'selected ';" << "\n";
        output << "            echo '>'.$Seed_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        $idxSeed = Explode( ':', $Seed_Select );" << "\n";
        output << "        $idxSeed = $idxSeed[0];" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='SVG' value='$SVG' />" << "\n";
        output << "            <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='MaxPPM' value='$MaxPPM' />" << "\n";
        output << "            <input type='hidden' name='SumPPM' value='$SumPPM' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--=================== Specific Length =====================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<input type=text name=Length size=3 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Length == '' )" << "\n";
        output << "        {" << "\n";
        output << "            echo 'Length';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $Length;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--==================== PPM Filtering ======================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Filter size=3 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $Filter == '' )" << "\n";
        output << "        {" << "\n";
        output << "            echo 'PPM';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $Filter;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "#<!--======================= Max PPM =========================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=MaxPPM size=4 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $MaxPPM == '' || $MaxPPM == 0 )" << "\n";
        output << "        {" << "\n";
        output << "            echo 'MaxPPM';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $MaxPPM;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        if( $MaxPPM == 'MaxPPM' || $MaxPPM == '' ) $MaxPPM = 0;" << "\n";
        output << "" << "\n";
        output << "#<!--======================= Sum PPM =========================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=SumPPM size=4 value=';" << "\n";
        output << "" << "\n";
        output << "        if( $SumPPM == '' || $SumPPM == 0 )" << "\n";
        output << "        {" << "\n";
        output << "            echo 'SumPPM';" << "\n";
        output << "        }" << "\n";
        output << "        else echo $SumPPM;" << "\n";
        output << "" << "\n";
        output << "        echo ' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        if( $SumPPM == 'SumPPM' || $SumPPM == '' ) $SumPPM = 0;" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Heatmap Color ======================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='SVG' value='$SVG' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Seed_Select' value='$Seed_Select' />" << "\n";
        output << "            <input type='hidden' name='Segment_Select' value='$Segment_Select' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--==================== SVG Converter ======================-->" << "\n";
        output << "" << "\n";
        output << "        if( $isRNAfold && $Seed_Select != '' && $Length != '' && $Length != 'Length' )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<select name=SVG onchange=this.form.submit();>';" << "\n";
        output << "            echo \"<option value='' \"   ; if($SVG=='')    echo 'selected'; echo '>WEB Mode</option>';" << "\n";
        output << "            echo \"<option value='Yes' \"; if($SVG=='Yes') echo 'selected'; echo '>SVG Mode</option>';" << "\n";
        output << "" << "\n";
        output << "            echo \"</select>" << "\n";
        output << "                <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "                <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "                <input type='hidden' name='MaxPPM' value='$MaxPPM' />" << "\n";
        output << "                <input type='hidden' name='SumPPM' value='$SumPPM' />" << "\n";
        output << "                <input type='hidden' name='Seed_Select' value='$Seed_Select' />" << "\n";
        output << "                <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "                <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "                <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "                <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "        }" << "\n";
        output << "        else $SVG = '';" << "\n";
        output << "" << "\n";
        output << "        $isSVG = $SVG == 'Yes' ? true : false ;" << "\n";
        output << "" << "\n";
        output << "        if( $isRNAfold && !$isSVG )" << "\n";
        output << "            echo 'MFE: '.$RNAfold[ 'MFE' ][0];" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Make TempFile ======================-->" << "\n";
        output << "" << "\n";
        output << "        $TSV_File_Name = Substr( $TSV_File, 0, ( Strlen( $TSV_File ) -4 ));" << "\n";
        output << "" << "\n";
        output << "        if( $Full_Sequc != '' && ( $Seg_End != 0 || $Data_Array[0][5] != '.' ))" << "\n";
        output << "        {" << "\n";
        output << "            $Temp = Tempnam( '/tmp', $TSV_File.'_'.$Annotation_Select.'_len'.$Length.'_ppm'.$Filter.'_JSON_' );" << "\n";
        output << "            $Ftemp = Fopen( $Temp, 'w' );" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp, \"{\\n\" );" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"Annotation\\\" : \\\"\".$Annotation_Select.\"\\\",\\n\" );" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"Sequence\\\"   : \\\"\"." << "\n";
        output << "                ( $Seg_End == 0 ? $Full_Sequc :" << "\n";
        output << "                    ( $Data_Array[0][5] == '.'" << "\n";
        output << "                        ? Substr( $Full_Sequc, $Seg_Start - " << extend_refseq << ", ( $Seg_End - $Seg_Start + " << ( extend_refseq * 2 ) << " ))" << "\n";
        output << "                        : Substr( $Full_Sequc, $Seg_Start, ( $Seg_End - $Seg_Start ))" << "\n";
        output << "                    )" << "\n";
        output << "                ).\"\\\",\\n\" );" << "\n";
        output << "" << "\n";
        output << "            if( $isRNAfold )" << "\n";
        output << "            {" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold0\\\"      : \\\"\".$RNAfold[ 'Fold0' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold1\\\"      : \\\"\".$RNAfold[ 'Fold1' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold2\\\"      : \\\"\".$RNAfold[ 'Fold2' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold3\\\"      : \\\"\".$RNAfold[ 'Fold3' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold4\\\"      : \\\"\".$RNAfold[ 'Fold4' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Fold5\\\"      : \\\"\".$RNAfold[ 'Fold5' ][0].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"MFE\\\"        : \".$RNAfold[ 'MFE' ][0].\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"MaxEntropy\\\" : \".$RNAfold[ 'MaxEntropy' ][0].\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"Entropy\\\"    : \\\"\".$RNAfold[ 'Entropy' ][0].\"\\\",\\n\" );" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"vEntropy\\\"   : [\" );" << "\n";
        output << "                for( $i = 0; $i < Count( $RNAfold[ 'vEntropy' ] ); ++$i )" << "\n";
        output << "                    Fwrite( $Ftemp, $RNAfold[ 'vEntropy' ][$i].( $i < Count( $RNAfold[ 'vEntropy' ] )-1 ? ',' : \"],\\n\" ));" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"xEntropy\\\"   : [\" );" << "\n";
        output << "                for( $i = 0; $i < Count( $RNAfold[ 'xEntropy' ] ); ++$i )" << "\n";
        output << "                    Fwrite( $Ftemp, $RNAfold[ 'xEntropy' ][$i].( $i < Count( $RNAfold[ 'xEntropy' ] )-1 ? ',' : \"],\\n\" ));" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"  \\\"yEntropy\\\"   : [\" );" << "\n";
        output << "                for( $i = 0; $i < Count( $RNAfold[ 'yEntropy' ] ); ++$i )" << "\n";
        output << "                    Fwrite( $Ftemp, $RNAfold[ 'yEntropy' ][$i].( $i < Count( $RNAfold[ 'yEntropy' ] )-1 ? ',' : \"],\\n\" ));" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"RefChr\\\"     : \\\"\".$Ref_Chr.\"\\\",\\n\" );" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"RefStart\\\"   : \".$Ref_Start.\",\\n\" );" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"RefStrand\\\"  : \\\"\".$Ref_Strand.\"\\\",\\n\" );" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp, \"  \\\"Reads\\\"      : [\\n\" );" << "\n";
        output << "" << "\n";
        output << "            $isFirst = true;" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Data_Array ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                if( $Seg_End != 0 && $Segm_Array[$i][0] != $Seg_Start && $Segm_Array[$i][1] != $Seg_End )" << "\n";
        output << "                    continue;" << "\n";
        output << "" << "\n";
        output << "                if( $Seed_Select != '' )" << "\n";
        output << "                {" << "\n";
        output << "                    $Seed = Explode( ':', $Seed_Select );" << "\n";
        output << "                    if( $Data_Array[$i][1] != $Seed[0] -2 ) continue;" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                if( !$isFirst ) Fwrite( $Ftemp, \",\\n\" );" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"    {\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"Index\\\"    : \"  .( $Data_Array[$i][1] - $Seg_Start + ( $Data_Array[0][5] == '.' ? " << extend_refseq << " : 0 )).  \",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"Length\\\"   : \"  .$Data_Array[$i][2].  \",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"PPM\\\"      : \"  .$Data_Array[$i][3].  \",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"isRMSK\\\"   : \\\"\".$Data_Array[$i][4].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"Arm\\\"      : \\\"\".$Data_Array[$i][5].\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"SegStart\\\" : \\\"\".  $Segm_Array[$i][0]  . \"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"SegEnd\\\"   : \\\"\".  $Segm_Array[$i][1]  . \"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"Tail\\\"     : \\\"\".( $Data_Array[$i][6] != \".\" ? $Data_Array[$i][6] : \"\" ).\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"MDtag\\\"    : \\\"\".( $Data_Array[$i][7] != \".\" ? $Data_Array[$i][7] : \"\" ).\"\\\",\\n\" );" << "\n";
        output << "                Fwrite( $Ftemp, \"      \\\"TCtag\\\"    : \\\"\".( $Data_Array[$i][8] != \".\" ? $Data_Array[$i][8] : \"\" ).\"\\\",\\n\" );" << "\n";
        output << "" << "\n";
        output << "                if( ( $Length != '' && $Length != 'Length' && $Length != $Data_Array[$i][2] ) ||" << "\n";
        output << "                    ( $Filter != '' && $Filter != 'PPM'    && $Filter >  $Data_Array[$i][3] ) )" << "\n";
        output << "                    Fwrite( $Ftemp, \"      \\\"Filter\\\"   : \\\"Y\\\"\\n\" );" << "\n";
        output << "                else" << "\n";
        output << "                    Fwrite( $Ftemp, \"      \\\"Filter\\\"   : \\\"N\\\"\\n\" );" << "\n";
        output << "" << "\n";
        output << "                Fwrite( $Ftemp, \"    }\" );" << "\n";
        output << "" << "\n";
        output << "                $isFirst = false;" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            Fwrite( $Ftemp, \"\\n  ]\\n\" );" << "\n";
        output << "            Fwrite( $Ftemp, \"}\\n\" );" << "\n";
        output << "            Fclose( $Ftemp );" << "\n";
        output << "" << "\n";
        output << "#<!--==================== Sq Align Plot ======================-->" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>" << "\n";
        output << "" << "\n";
        output << "                var color_map = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);" << "\n";
        output << "" << "\n";
        output << "                var datas = Array();" << "\n";
        output << "                var exprs = Array();" << "\n";
        output << "                var etpys = Array();" << "\n";
        output << "                var seeds = Array();" << "\n";
        output << "" << "\n";
        output << "                var anno = '';" << "\n";
        output << "                var sequ = '';" << "\n";
        output << "" << "\n";
        output << "                var colorA = 'green';" << "\n";
        output << "                var colorC = 'blue';" << "\n";
        output << "                var colorG = 'goldenrod';" << "\n";
        output << "                var colorU = 'red';" << "\n";
        output << "                var colorP = 'darkgray';\";" << "\n";
        output << "" << "\n";
        output << "            if( $isRNAfold )" << "\n";
        output << "            {" << "\n";
        output << "                echo \"" << "\n";
        output << "                var shift_top2 = 38;" << "\n";
        output << "                var fold_height = 200;\";" << "\n";
        output << "            }" << "\n";
        output << "            else" << "\n";
        output << "            {" << "\n";
        output << "                echo \"" << "\n";
        output << "                var shift_top2 = 0;" << "\n";
        output << "                var fold_height = 0;\";" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo \"" << "\n";
        output << "                var exp_height = 120;" << "\n";
        output << "                var seq_height = 50;" << "\n";
        output << "                var read_height = 20;" << "\n";
        output << "                var spc_num = 18;" << "\n";
        output << "                var expr_max = 0;" << "\n";
        output << "                var shift_top1 = 30 + fold_height;" << "\n";
        output << "" << "\n";
        output << "                var seg_start = $Seg_Start;" << "\n";
        output << "                var seg_end   = $Seg_End;" << "\n";
        output << "" << "\n";
        output << "                $.getJSON( '$Temp', function( json )" << "\n";
        output << "                {" << "\n";
        output << "                    anno = json[ 'Annotation' ];" << "\n";
        output << "                    sequ = json[ 'Sequence'   ];" << "\n";
        output << "" << "\n";
        output << "                    for( i in json[ 'vEntropy' ] )" << "\n";
        output << "                        etpys[i] = json[ 'vEntropy' ][i];" << "\n";
        output << "" << "\n";
        output << "                    var expr_array = Array();" << "\n";
        output << "                    var heat_array = Array();" << "\n";
        output << "" << "\n";
        output << "                    var heat_max = 0;" << "\n";
        output << "                    var lb_width = 1;" << "\n";
        output << "" << "\n";
        output << "                    var start = 0;" << "\n";
        output << "                    var end   = 0;" << "\n";
        output << "                    var extra = 0;" << "\n";
        output << "" << "\n";
        output << "                    $.each( json[ 'Reads' ], function( idx, read )" << "\n";
        output << "                    {" << "\n";
        output << "                        end = read[ 'Index' ] + read[ 'Length' ] + read[ 'Tail' ].length;" << "\n";
        output << "" << "\n";
        output << "                        if(( end - json[ 'Sequence' ].length ) > 0 &&" << "\n";
        output << "                           ( end - json[ 'Sequence' ].length ) > extra )" << "\n";
        output << "                        {" << "\n";
        output << "                            extra = end - json[ 'Sequence' ].length;" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        if( !(( read[ 'Index' ] + 4 ) in seeds ))" << "\n";
        output << "                        {" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ] = Array();" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][ 'Name' ]" << "\n";
        output << "                                = anno + \".( $Data_Array[0][5] != '.' ? ( \"'-' + read[ 'Arm' ]\" ) : \"''\" ).\" + '_' + sequ.substr(( read[ 'Index' ] + 1 ), 7 );" << "\n";
        output << "" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][ 'GMPM' ] = 0; " << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][  'GM'  ] = 0;" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][  'PM'  ] = 0;" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][ 'Counts' ] = 0;" << "\n";
        output << "                            seeds[ read[ 'Index' ] + 4 ][ 'Length' ] = 0;" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        seeds[ read[ 'Index' ] + 4 ][ 'GMPM' ] += read[ 'PPM' ];" << "\n";
        output << "                        seeds[ read[ 'Index' ] + 4 ][  'GM'  ] += ( read[ 'Tail' ] == '' ? read[ 'PPM' ] : 0 );" << "\n";
        output << "                        seeds[ read[ 'Index' ] + 4 ][  'PM'  ] += ( read[ 'Tail' ] != '' ? read[ 'PPM' ] : 0 );" << "\n";
        output << "                        seeds[ read[ 'Index' ] + 4 ][ 'Counts' ]++;" << "\n";
        output << "                        seeds[ read[ 'Index' ] + 4 ][ 'Length' ]" << "\n";
        output << "                            = (( read[ 'Length' ] + read[ 'Tail' ].length ) > seeds[ read[ 'Index' ] + 4 ][ 'Length' ]" << "\n";
        output << "                            ?  ( read[ 'Length' ] + read[ 'Tail' ].length )" << "\n";
        output << "                            :    seeds[ read[ 'Index' ] + 4 ][ 'Length' ] )" << "\n";
        output << "                            ;" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Sequence' ].length + extra; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        expr_array[ i ] = Array();" << "\n";
        output << "                        expr_array[ i ][ 'PPM' ] = 0;" << "\n";
        output << "                        expr_array[ i ][  'P'  ] = 0;" << "\n";
        output << "                        expr_array[ i ][  'A'  ] = 0;" << "\n";
        output << "                        expr_array[ i ][  'C'  ] = 0;" << "\n";
        output << "                        expr_array[ i ][  'G'  ] = 0;" << "\n";
        output << "                        expr_array[ i ][  'U'  ] = 0;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $.each( json[ 'Reads' ], function( idx, read )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( read[ 'Filter' ] == 'Y' ) return;" << "\n";
        output << "" << "\n";
        output << "                        datas[ idx ] = read;" << "\n";
        output << "                        heat_array[ idx ] = read[ 'PPM' ];" << "\n";
        output << "" << "\n";
        output << "                        start = read[ 'Index' ];" << "\n";
        output << "                        end   = read[ 'Index' ] + read[ 'Length' ];" << "\n";
        output << "                        pos   = '';" << "\n";
        output << "" << "\n";
        output << "                        for( var i = start; i < ( end + read[ 'Tail' ].length ); ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            expr_array[ i ][ 'PPM' ] += read[ 'PPM' ];" << "\n";
        output << "" << "\n";
        output << "                            if( i < end ) expr_array[ i ][  'P'  ] += read[ 'PPM' ];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'MDtag' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( !$.isNumeric( read[ 'MDtag' ][i] ))" << "\n";
        output << "                            {" << "\n";
        output << "                                expr_array[ start + parseInt( pos )][ read[ 'MDtag' ][i] ] += read[ 'PPM' ];" << "\n";
        output << "                                pos = '';" << "\n";
        output << "                            }" << "\n";
        output << "                            else pos += read[ 'MDtag' ][i];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'TCtag' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( !$.isNumeric( read[ 'TCtag' ][i] ))" << "\n";
        output << "                            {" << "\n";
        output << "                                expr_array[ start + parseInt( pos )][ read[ 'TCtag' ][i] ] += read[ 'PPM' ];" << "\n";
        output << "                                pos = '';" << "\n";
        output << "                            }" << "\n";
        output << "                            else pos += read[ 'TCtag' ][i];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'Tail' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'A' ) expr_array[ i + end ][ 'A' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'C' ) expr_array[ i + end ][ 'C' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'G' ) expr_array[ i + end ][ 'G' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'U' ) expr_array[ i + end ][ 'U' ] += read[ 'PPM' ];" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( expr_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( expr_max < expr_array[ idx ][ 'PPM' ] ) expr_max = expr_array[ idx ][ 'PPM' ];" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    if( $SumPPM != 0 && $SumPPM != 'SumPPM' )" << "\n";
        output << "                        expr_max = $SumPPM;" << "\n";
        output << "" << "\n";
        output << "                    $.each( expr_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        expr_array[ idx ][ 'PPM' ] = expr_array[ idx ][ 'PPM' ] / expr_max * 100;" << "\n";
        output << "" << "\n";
        output << "                        var total = expr_array[ idx ][ 'P' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'A' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'C' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'G' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'U' ];" << "\n";
        output << "" << "\n";
        output << "                        if( total != 0 )" << "\n";
        output << "                        {" << "\n";
        output << "                            expr_array[ idx ][ 'P' ] = expr_array[ idx ][ 'P' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'A' ] = expr_array[ idx ][ 'A' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'C' ] = expr_array[ idx ][ 'C' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'G' ] = expr_array[ idx ][ 'G' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'U' ] = expr_array[ idx ][ 'U' ] / total;" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( heat_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( heat_max < heat_array[ idx ] )" << "\n";
        output << "                            heat_max = heat_array[ idx ];" << "\n";
        output << "                    });\";" << "\n";
        output << "" << "\n";
        output << "            if( $MaxPPM != 0 && $MaxPPM != 'MaxPPM' )" << "\n";
        output << "                echo \"                    heat_max = $MaxPPM;\";" << "\n";
        output << "" << "\n";
        output << "            echo \"" << "\n";
        output << "                    $.each( heat_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        heat_array[ idx ] = heat_array[ idx ] / heat_max * 100;" << "\n";
        output << "                        if( heat_array[ idx ] == 0 ) heat_array[ idx ] = 0.01;" << "\n";
        output << "                    });\";" << "\n";
        output << "" << "\n";
        output << "            if( $isSVG ) echo \"" << "\n";
        output << "                    spc_num = 9;" << "\n";
        output << "                    fold_height = 100;" << "\n";
        output << "                    read_height = 12;" << "\n";
        output << "" << "\n";
        output << "                    $( 'body' ).append( \\\"<svg id='svg'></svg>\\\" );" << "\n";
        output << "                    var svg = d3.select( '#svg' )" << "\n";
        output << "                        .attr( 'xmlns', 'http://www.w3.org/2000/svg' )" << "\n";
        output << "                        .attr( 'version', '1.1' )" << "\n";
        output << "                        .attr( 'width', (( json[ 'Sequence' ].length + extra + lb_width + 7 )* spc_num ) + 'px' )" << "\n";
        output << "                        ;" << "\n";
        output << "            \";" << "\n";
        output << "" << "\n";
        output << "            else echo \"" << "\n";
        output << "                    $( 'body' ).append( \\\"<div id='chart'></div>\\\" );" << "\n";
        output << "                    $( '#chart' ).css({" << "\n";
        output << "                        'margin-top': '10px'," << "\n";
        output << "                        });\";" << "\n";
        output << "" << "\n";
        output << "            if( $isRNAfold )" << "\n";
        output << "            {" << "\n";
        output << "                if( $isSVG ) echo \"" << "\n";
        output << "                    var foldchart = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'foldchart' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "                \";" << "\n";
        output << "                else echo\"" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='foldchart'></div>\\\" );" << "\n";
        output << "                    $( '#foldchart' ).css({" << "\n";
        output << "                        'height': fold_height," << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        });" << "\n";
        output << "                \";" << "\n";
        output << "" << "\n";
        output << "                echo \"" << "\n";
        output << "                    var foldidx = 0;" << "\n";
        output << "                    var foldidxs = [];" << "\n";
        output << "" << "\n";
        output << "                    for( var f = 1; f <= 5; f++ )" << "\n";
        output << "                        foldidxs[f] = [];" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Fold1' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        for( var f = 1; f <= 2; f++ )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( json[ 'Fold' + f ].charAt(i) != ' ' )" << "\n";
        output << "                            {" << "\n";
        output << "                                foldidx++;" << "\n";
        output << "                                foldidxs[f][i] = foldidx;" << "\n";
        output << "                            }" << "\n";
        output << "                            else foldidxs[f][i] = 0;" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Fold3' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( json[ 'Fold3' ].charAt(i) != ' ' && json[ 'Fold3' ].charAt(i) != '|' && json[ 'Fold3' ].charAt(i) != ':' && json[ 'Fold3' ].charAt(i) != ')' )" << "\n";
        output << "                        {" << "\n";
        output << "                            foldidx++;" << "\n";
        output << "                            foldidxs[3][i] = foldidx;" << "\n";
        output << "                        }" << "\n";
        output << "                        else foldidxs[3][i] = 0;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    for( var i = json[ 'Fold5' ].length -1; i >= 0; --i )" << "\n";
        output << "                    {" << "\n";
        output << "                        for( var f = 4; f <= 5; f++ )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( json[ 'Fold' + f ].charAt(i) != ' ' )" << "\n";
        output << "                            {" << "\n";
        output << "                                foldidx++;" << "\n";
        output << "                                foldidxs[f][i] = foldidx;" << "\n";
        output << "                            }" << "\n";
        output << "                            else foldidxs[f][i] = 0;" << "\n";
        output << "                        }" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    var right_shift = 0;" << "\n";
        output << "                    var left_shift = (( json[ 'Sequence' ].length * spc_num ) /2 ) - (( json[ 'Fold3' ].length + 18 ) * spc_num / 2 );" << "\n";
        output << "" << "\n";
        output << "                    for( var f = 1; f <= 5; f++ )" << "\n";
        output << "                    {" << "\n";
        output << "                        var foldx = 'Fold' + f;" << "\n";
        output << "                        var shift = (( json[ 'Sequence' ].length * spc_num ) /2 ) - (( json[ foldx ].length + 18 ) * spc_num / 2 );" << "\n";
        output << "                \";" << "\n";
        output << "" << "\n";
        output << "                if( $isSVG ) echo \"" << "\n";
        output << "                        var foldg = foldchart.append( 'g' )" << "\n";
        output << "                            .attr( 'id', foldx )" << "\n";
        output << "                            .attr( 'transform', 'translate(' + ( shift + 40 ) + ',' + ( f * 18 ) + ')' )" << "\n";
        output << "                            ;" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < json[ foldx ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( json[ foldx ].charAt(i) == ' ' )" << "\n";
        output << "                                continue;" << "\n";
        output << "" << "\n";
        output << "                            y = 0;" << "\n";
        output << "                            x = ( i * spc_num );" << "\n";
        output << "" << "\n";
        output << "                            if( f == 3 )" << "\n";
        output << "                            {" << "\n";
        output << "                                y = -2;" << "\n";
        output << "                                x = ( i * spc_num ) + 2.5;" << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            if( foldidxs[f][i] >= ( $idxSeed - 1 ) && foldidxs[f][i] < ( $idxSeed + $Length - 1 ))" << "\n";
        output << "                                foldg.append( 'rect' )" << "\n";
        output << "                                    .attr( 'x', x-0.5 )" << "\n";
        output << "                                    .attr( 'y', -12 )" << "\n";
        output << "                                    .attr( 'width', '10px' )" << "\n";
        output << "                                    .attr( 'height', '14px' )" << "\n";
        output << "                                    .attr( 'fill', 'gold' )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                            foldg.append( 'text' )" << "\n";
        output << "                                .attr( 'x', x )" << "\n";
        output << "                                .attr( 'y', y )" << "\n";
        output << "                                .text( json[ foldx ].charAt(i) )" << "\n";
        output << "                                ;" << "\n";
        output << "                        }" << "\n";
        output << "                \";" << "\n";
        output << "                else echo\"" << "\n";
        output << "                        $( '#foldchart' ).append( \\\"<div id='\\\" + foldx + \\\"'></div>\\\" );" << "\n";
        output << "                        $( '#' + foldx ).css({" << "\n";
        output << "                            'position': 'relative'," << "\n";
        output << "                            'left': shift + 'px'," << "\n";
        output << "                            'width': ( json[ 'Sequence' ].length * spc_num ) - shift + 'px'," << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '30px'," << "\n";
        output << "                            });" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < json[ foldx ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            $( '#' + foldx ).append( \\\"<div id='fold\\\" + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] ) + \\\"' title='fold\\\" + foldidxs[f][i] + \\\"'>\\\" + json[ foldx ].charAt(i) + '</div>' );" << "\n";
        output << "                            $( '#fold' + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] )).css({ 'white-space': 'pre' });" << "\n";
        output << "                            if( i != json[ foldx ].length-1 )" << "\n";
        output << "                                $( '#fold' + ( foldidxs[f][i] == 0 ? ( f*1000+i ) : foldidxs[f][i] )).css({ 'float': 'left' });" << "\n";
        output << "                        }\";" << "\n";
        output << "" << "\n";
        output << "                echo '                    }';" << "\n";
        output << "" << "\n";
        output << "                if( $isSVG ) echo \"" << "\n";
        output << "                    var info = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'info' )" << "\n";
        output << "                        .attr( 'transform', 'translate(' + ( left_shift + json[ 'Fold3' ].length * spc_num + 80 ) + ', 26)' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    info.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', -20 )" << "\n";
        output << "                        .attr( 'x2', -20 )" << "\n";
        output << "                        .attr( 'y1', -20 )" << "\n";
        output << "                        .attr( 'y2',  70 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    info.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 0 )" << "\n";
        output << "                        .attr( 'y', 0 )" << "\n";
        output << "                        .text( 'miRNA: \".$Annotation_Select.\"' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    info.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 12.6 )" << "\n";
        output << "                        .attr( 'y', 18 )" << "\n";
        output << "                        .text( 'MFE: \".$RNAfold[ 'MFE' ][0].\"' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    info.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 9.9 )" << "\n";
        output << "                        .attr( 'y', 36 )" << "\n";
        output << "                        .text( 'Seed: \".Substr( $Seed_Select, -7 ).\"' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    info.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 1.3 )" << "\n";
        output << "                        .attr( 'y', 54 )" << "\n";
        output << "                        .text( 'Length: \".$Length.\" nt' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    svg.append( 'rect' )" << "\n";
        output << "                        .attr( 'id', 'selectbox' )" << "\n";
        output << "                        .attr( 'height', '78px' )" << "\n";
        output << "                        .attr( 'x', ( $idxSeed -2 ) * spc_num + 26 )" << "\n";
        output << "                        .attr( 'y', fold_height + 8 )" << "\n";
        output << "                        .attr( 'width', ( $Length * spc_num ) + 'px' )" << "\n";
        output << "                        .attr( 'fill', 'gold' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    var entropy = svg.append( 'g' )" << "\n";
        output << "                       .attr( 'id', 'entropy' )" << "\n";
        output << "                       .attr( 'transform', 'translate( 20,' + ( fold_height + 10 ) + ')' )" << "\n";
        output << "                       .attr( 'font-family', 'Arial' )" << "\n";
        output << "                       .attr( 'font-size', '12px' )" << "\n";
        output << "                       ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'text' )" << "\n";
        output << "                        .attr( 'x', -42 )" << "\n";
        output << "                        .attr( 'y', -1 )" << "\n";
        output << "                        .attr( 'transform', 'rotate(270)' )" << "\n";
        output << "                        .text( 'Entropy' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 5 )" << "\n";
        output << "                        .attr( 'x2', 5 )" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 46 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'polyline' )" << "\n";
        output << "                        .attr( 'points', json[ 'Entropy' ])" << "\n";
        output << "                        .attr( 'fill', 'none' )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', 2 )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', (( json[ 'Sequence' ].length + extra ) * spc_num ))" << "\n";
        output << "                        .attr( 'x2', (( json[ 'Sequence' ].length + extra ) * spc_num ))" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 46 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', (( json[ 'Sequence' ].length + extra ) * spc_num ))" << "\n";
        output << "                        .attr( 'x2', (( json[ 'Sequence' ].length + extra ) * spc_num + 5 ))" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 0 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'text' )" << "\n";
        output << "                        .attr( 'x', (( json[ 'Sequence' ].length + extra ) * spc_num + 7 ))" << "\n";
        output << "                        .attr( 'y', 4 )" << "\n";
        output << "                        .text( json[ 'MaxEntropy' ].toFixed(2) )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', (( json[ 'Sequence' ].length + extra ) * spc_num ))" << "\n";
        output << "                        .attr( 'x2', (( json[ 'Sequence' ].length + extra ) * spc_num + 5 ))" << "\n";
        output << "                        .attr( 'y1', 46 )" << "\n";
        output << "                        .attr( 'y2', 46 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    entropy.append( 'text' )" << "\n";
        output << "                        .attr( 'x', (( json[ 'Sequence' ].length + extra ) * spc_num + 7 ))" << "\n";
        output << "                        .attr( 'y', 50 )" << "\n";
        output << "                        .text( '0.00' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                     var squence = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'squence' )" << "\n";
        output << "                        .attr( 'transform', 'translate( 26,' + ( fold_height + 68 ) + ')' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    var seqs = squence.append( 'g' )" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Sequence' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        seqs.append( 'text' )" << "\n";
        output << "                            .attr( 'y', 2 )" << "\n";
        output << "                            .attr( 'x', ( i * spc_num ))" << "\n";
        output << "                            .text( json[ 'Sequence' ].charAt(i) )" << "\n";
        output << "                            ;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    var fold = squence.append( 'g' )" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Fold0' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        y = 15;" << "\n";
        output << "" << "\n";
        output << "                        if( json[ 'Fold0' ].charAt(i) == '.' )" << "\n";
        output << "                            y = 12;" << "\n";
        output << "" << "\n";
        output << "                        fold.append( 'text' )" << "\n";
        output << "                            .attr( 'y', y )" << "\n";
        output << "                            .attr( 'x', ( i * spc_num + 3 ))" << "\n";
        output << "                            .text( json[ 'Fold0' ].charAt(i) )" << "\n";
        output << "                            ;" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    var first_idx = -1;" << "\n";
        output << "                    var last_idx = 0;" << "\n";
        output << "" << "\n";
        output << "                    var expression = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'expression' )" << "\n";
        output << "                        .attr( 'transform', 'translate( 26,' + ( fold_height + 100 ) + ')' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    $.each( expr_array, function( idx, expr )" << "\n";
        output << "                    {" << "\n";
        output << "                        exprs[ idx ] = expr;" << "\n";
        output << "" << "\n";
        output << "                        if( expr[ 'U' ] + expr[ 'G' ] + expr[ 'C' ] + expr[ 'A' ] + expr[ 'P' ] != 0 )" << "\n";
        output << "                        {" << "\n";
        output << "                            last_idx = idx;" << "\n";
        output << "                            if( first_idx == -1 )" << "\n";
        output << "                                first_idx = idx;" << "\n";
        output << "" << "\n";
        output << "                            var expressiong = expression.append( 'g' )" << "\n";
        output << "                                .attr( 'transform', 'translate(' + ( idx * spc_num ) + ',' + (( 100 - ( expr[ 'PPM' ] * ( expr[ 'U' ] + expr[ 'G' ] + expr[ 'C' ] + expr[ 'A' ] + expr[ 'P' ] )) ) /2  ) + ')' )" << "\n";
        output << "                                ;" << "\n";
        output << "" << "\n";
        output << "                            if( expr[ 'U' ] != 0 )" << "\n";
        output << "                                expressiong.append( 'rect' )" << "\n";
        output << "                                    .attr( 'id', 'U' )" << "\n";
        output << "                                    .attr( 'x', 0 )" << "\n";
        output << "                                    .attr( 'y', 0 )" << "\n";
        output << "                                    .attr( 'height', (( expr[ 'PPM' ] * expr[ 'U' ]) /2 ) + 'px' )" << "\n";
        output << "                                    .attr( 'width', spc_num + 'px' )" << "\n";
        output << "                                    .attr( 'fill', colorU )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                            if( expr[ 'G' ] != 0 )" << "\n";
        output << "                                expressiong.append( 'rect' )" << "\n";
        output << "                                    .attr( 'id', 'G' )" << "\n";
        output << "                                    .attr( 'x', 0 )" << "\n";
        output << "                                    .attr( 'y', ( expr[ 'PPM' ] * expr[ 'U' ]) /2 )" << "\n";
        output << "                                    .attr( 'height', (( expr[ 'PPM' ] * expr[ 'G' ]) /2 ) + 'px' )" << "\n";
        output << "                                    .attr( 'width', spc_num + 'px' )" << "\n";
        output << "                                    .attr( 'fill', colorG )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                            if( expr[ 'C' ] != 0 )" << "\n";
        output << "                                expressiong.append( 'rect' )" << "\n";
        output << "                                    .attr( 'id', 'C' )" << "\n";
        output << "                                    .attr( 'x', 0 )" << "\n";
        output << "                                    .attr( 'y', ( expr[ 'PPM' ] * ( expr[ 'U' ] + expr[ 'G' ])) /2 )" << "\n";
        output << "                                    .attr( 'height', (( expr[ 'PPM' ] * expr[ 'C' ]) /2 ) + 'px' )" << "\n";
        output << "                                    .attr( 'width', spc_num + 'px' )" << "\n";
        output << "                                    .attr( 'fill', colorC )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                            if( expr[ 'A' ] != 0 )" << "\n";
        output << "                                expressiong.append( 'rect' )" << "\n";
        output << "                                    .attr( 'id', 'A' )" << "\n";
        output << "                                    .attr( 'x', 0 )" << "\n";
        output << "                                    .attr( 'y', ( expr[ 'PPM' ] * ( expr[ 'U' ] + expr[ 'G' ] + expr[ 'C' ])) /2 )" << "\n";
        output << "                                    .attr( 'height', (( expr[ 'PPM' ] * expr[ 'A' ]) /2 ) + 'px' )" << "\n";
        output << "                                    .attr( 'width', spc_num + 'px' )" << "\n";
        output << "                                    .attr( 'fill', colorA )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                            if( expr[ 'P' ] != 0 )" << "\n";
        output << "                                expressiong.append( 'rect' )" << "\n";
        output << "                                    .attr( 'id', 'P' )" << "\n";
        output << "                                    .attr( 'x', 0 )" << "\n";
        output << "                                    .attr( 'y', ( expr[ 'PPM' ] * ( expr[ 'U' ] + expr[ 'G' ] + expr[ 'C' ] + expr[ 'A' ])) /2 )" << "\n";
        output << "                                    .attr( 'height', (( expr[ 'PPM' ] * expr[ 'P' ]) /2 ) + 'px' )" << "\n";
        output << "                                    .attr( 'width', spc_num + 'px' )" << "\n";
        output << "                                    .attr( 'fill', colorP )" << "\n";
        output << "                                    ;" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    var explegend = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'expression_legend' )" << "\n";
        output << "                        .attr( 'transform', 'translate(' + ((( last_idx +1 ) * spc_num ) + 29 ) + ',' + ( fold_height + 100 ) + ')' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    explegend.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 0 )" << "\n";
        output << "                        .attr( 'x2', 0 )" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 50 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    explegend.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 0 )" << "\n";
        output << "                        .attr( 'x2', 5 )" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 0 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    explegend.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 7 )" << "\n";
        output << "                        .attr( 'y', 4 )" << "\n";
        output << "                        .text( expr_max.toFixed(0) )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    explegend.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 0 )" << "\n";
        output << "                        .attr( 'x2', 5 )" << "\n";
        output << "                        .attr( 'y1', 50 )" << "\n";
        output << "                        .attr( 'y2', 50 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    explegend.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 7 )" << "\n";
        output << "                        .attr( 'y', 54 )" << "\n";
        output << "                        .text( '0' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    var title = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'title' )" << "\n";
        output << "                        .attr( 'transform', 'translate(' + ( first_idx * spc_num + 26 ) + ',' + ( fold_height + 170 ) + ')' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    title.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 0 )" << "\n";
        output << "                        .attr( 'x2', (( last_idx - first_idx + 1 - \".Strlen( Substr( $TSV_File, 0, -4 )).\" ) /2 * spc_num -5 ))" << "\n";
        output << "                        .attr( 'y1', -4 )" << "\n";
        output << "                        .attr( 'y2', -4 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    title.append( 'text' )" << "\n";
        output << "                        .attr( 'x', (( last_idx - first_idx + 1.25 - \".Strlen( Substr( $TSV_File, 0, -4 )).\" ) /2 * spc_num ))" << "\n";
        output << "                        .attr( 'y', 0 )" << "\n";
        output << "                        .text( '\".Substr( $TSV_File, 0, -4 ).\"' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    title.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', ( ( last_idx - first_idx + 1 + \".Strlen( Substr( $TSV_File, 0, -4 )).\" ) /2 * spc_num +5 ))" << "\n";
        output << "                        .attr( 'x2', (( last_idx - first_idx + 1 ) * spc_num ))" << "\n";
        output << "                        .attr( 'y1', -4 )" << "\n";
        output << "                        .attr( 'y2', -4 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    var legend = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'legend' )" << "\n";
        output << "                        .attr( 'transform', 'translate(' + ((( last_idx +1 ) * spc_num ) + 29 ) + ',' + ( fold_height + 180 ) + ')' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    var defs = legend.append( 'defs' )" << "\n";
        output << "                        .append( 'linearGradient' )" << "\n";
        output << "                        .attr( 'id', 'gradient' )" << "\n";
        output << "                        .attr( 'x1', '100%' )" << "\n";
        output << "                        .attr( 'y1', '100%' )" << "\n";
        output << "                        .attr( 'x2', '100%' )" << "\n";
        output << "                        .attr( 'y2', '0%' )" << "\n";
        output << "                        .attr( 'spreadMethod', 'pad' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    defs.append( 'stop' )" << "\n";
        output << "                        .attr( 'offset', '0%' )" << "\n";
        output << "                        .attr( 'stop-color', color_map( 0 ))" << "\n";
        output << "                        .attr( 'stop-opacity', 1 )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    defs.append( 'stop' )" << "\n";
        output << "                        .attr( 'offset', '100%' )" << "\n";
        output << "                        .attr( 'stop-color', color_map( 100 ))" << "\n";
        output << "                        .attr( 'stop-opacity', 1 )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    legend.append( 'rect' )" << "\n";
        output << "                        .attr( 'x', 0 )" << "\n";
        output << "                        .attr( 'y', 0 )" << "\n";
        output << "                        .attr( 'width', 8 )" << "\n";
        output << "                        .attr( 'height', 70 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        .attr( 'fill', 'url( #gradient )' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    legend.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 8 )" << "\n";
        output << "                        .attr( 'x2', 12 )" << "\n";
        output << "                        .attr( 'y1', 0 )" << "\n";
        output << "                        .attr( 'y2', 0 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    legend.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 14 )" << "\n";
        output << "                        .attr( 'y', 5 )" << "\n";
        output << "                        .text( heat_max.toFixed(0) )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    legend.append( 'line' )" << "\n";
        output << "                        .attr( 'x1', 8 )" << "\n";
        output << "                        .attr( 'x2', 12 )" << "\n";
        output << "                        .attr( 'y1', 70 )" << "\n";
        output << "                        .attr( 'y2', 70 )" << "\n";
        output << "                        .attr( 'stroke', 'black' )" << "\n";
        output << "                        .attr( 'stroke-width', '1px' )" << "\n";
        output << "                        ;" << "\n";
        output << "" << "\n";
        output << "                    legend.append( 'text' )" << "\n";
        output << "                        .attr( 'x', 14 )" << "\n";
        output << "                        .attr( 'y', 75 )" << "\n";
        output << "                        .text( '0' )" << "\n";
        output << "                        ;" << "\n";
        output << "                \";" << "\n";
        output << "            }" << "\n";
        output << "            else if( $isSVG )" << "\n";
        output << "            {" << "\n";
        output << "                // svg for no rnaFold" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            if( !$isSVG )" << "\n";
        output << "            {" << "\n";
        output << "                echo \"" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='expression'></div>\\\" );" << "\n";
        output << "                    $( '#expression' ).css({" << "\n";
        output << "                        'height': '120px'," << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#expression' ).append( \\\"<div id='exprchart'></div>\\\" );" << "\n";
        output << "                    $( '#exprchart' ).css({" << "\n";
        output << "                        'height': '104px'," << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        });\";" << "\n";
        output << "" << "\n";
        output << "                if( $isRNAfold ) echo \"" << "\n";
        output << "                    $( '#expression' ).append( $(\\\"<svg id='etpchart' xmlns='http://www.w3.org/2000/svg'><polyline id='entropy' points='\\\" + json[ 'Entropy' ] + \\\"'/></svg>\\\"))" << "\n";
        output << "                    $( '#etpchart' ).css({" << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'z-index': '-1'," << "\n";
        output << "                        'height': '100px'," << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        'top': 7 + shift_top1 + 'px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#entropy' ).css({" << "\n";
        output << "                        'fill': 'none'," << "\n";
        output << "                        'stroke': 'Paleturquoise'," << "\n";
        output << "                        'stroke-width': 2" << "\n";
        output << "                        });\";" << "\n";
        output << "" << "\n";
        output << "                echo \"" << "\n";
        output << "                    $.each( expr_array, function( idx, expr )" << "\n";
        output << "                    {" << "\n";
        output << "                        exprs[ idx ] = expr;" << "\n";
        output << "" << "\n";
        output << "                        $( '#exprchart' ).append( \\\"<div id='exprtitle\\\" + idx + \\\"' title='expr\\\" + idx + \\\"'></div>\\\" );" << "\n";
        output << "                        $( '#exprtitle' + idx ).css({" << "\n";
        output << "                            'height': '100px'," << "\n";
        output << "                            'width': spc_num + 'px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "" << "\n";
        output << "                        $( '#exprtitle' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).css({" << "\n";
        output << "                            'margin-top': ( 100 - expr[ 'PPM' ] ) + 'px'," << "\n";
        output << "                            'width': spc_num + 'px'," << "\n";
        output << "                            'position': 'relative'," << "\n";
        output << "                            'float': 'left'," << "\n";
        output << "                            'z-index': '-1'" << "\n";
        output << "                            });" << "\n";
        output << "" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"U'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"G'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"C'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"A'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"P'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                        $( '#expr' + idx + 'U' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'U' ]) + 'px', 'background': colorU });" << "\n";
        output << "                        $( '#expr' + idx + 'G' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'G' ]) + 'px', 'background': colorG });" << "\n";
        output << "                        $( '#expr' + idx + 'C' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'C' ]) + 'px', 'background': colorC });" << "\n";
        output << "                        $( '#expr' + idx + 'A' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'A' ]) + 'px', 'background': colorA });" << "\n";
        output << "                        $( '#expr' + idx + 'P' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'P' ]) + 'px', 'background': colorP });" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $( '#exprchart' ).append( \\\"<div id='exprlabel'></div>\\\" );" << "\n";
        output << "                    $( '#exprlabel' ).css({" << "\n";
        output << "                        'width': '1px'," << "\n";
        output << "                        'height': '100px'," << "\n";
        output << "                        'background': 'black'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });\";" << "\n";
        output << "" << "\n";
        output << "                if( $isRNAfold ) echo \"" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabeltop'>-\\\" + json[ 'MaxEntropy' ].toFixed(2) + '</div>' );\";" << "\n";
        output << "                else echo \"" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabeltop'>-\\\" + expr_max.toFixed(2) + '</div>' );\";" << "\n";
        output << "" << "\n";
        output << "                echo \"" << "\n";
        output << "                    $( '#exprlabeltop' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': 0 + shift_top1 + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabeldown'>-0.00</div>\\\" );" << "\n";
        output << "                    $( '#exprlabeldown' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': 100 + shift_top1 + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#expression' ).append( \\\"<div id='basecount1'></div>\\\" );" << "\n";
        output << "                    $( '#basecount1' ).css({ 'height': '10px' });" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 1; i < json[ 'Sequence' ].length + extra; ++i ) if( i % 5 == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#basecount1' ).append( \\\"<div id='base1-\\\" + i + \\\"' align='right'>\\\" + ( i == 5 ? '05' : i ) + '</div>' );" << "\n";
        output << "                        $( '#base1-' + i ).css({" << "\n";
        output << "                            'width': ( 5 * spc_num ) + 'px'," << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '10px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='sequence'></div>\\\" );" << "\n";
        output << "                    $( '#sequence' ).css({" << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        'height': seq_height + 'px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#sequence' ).append( \\\"<div id='selectbox'></div>\\\" );" << "\n";
        output << "                    $( '#selectbox' ).css({" << "\n";
        output << "                        'border': '5px solid gold'," << "\n";
        output << "                        'border-radius': '12px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'pointer-events': 'none',\";" << "\n";
        output << "" << "\n";
        output << "                if( $isRNAfold ) echo \"" << "\n";
        output << "                        'height': '68px',\";" << "\n";
        output << "                else echo \"" << "\n";
        output << "                        'height': '28px',\";" << "\n";
        output << "" << "\n";
        output << "                echo \"" << "\n";
        output << "                        'display': 'none'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#sequence' ).append( \\\"<div id='seqs'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                    $( '#seqs' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '30px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Sequence' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( i in seeds )" << "\n";
        output << "                        {" << "\n";
        output << "                            var MDarray = Array();" << "\n";
        output << "                            var MDseed = '';" << "\n";
        output << "                            var RMSK = '';" << "\n";
        output << "" << "\n";
        output << "                            $.each( json[ 'Reads' ], function( idx, read )" << "\n";
        output << "                            {" << "\n";
        output << "                                if( read[ 'Index' ] != i ) return;" << "\n";
        output << "                                RMSK = read[ 'isRMSK' ];" << "\n";
        output << "                                var pos = '';" << "\n";
        output << "" << "\n";
        output << "                                for( var j = 0; j < read[ 'MDtag' ].length; ++j )" << "\n";
        output << "                                {" << "\n";
        output << "                                    if( !$.isNumeric( read[ 'MDtag' ][j] ))" << "\n";
        output << "                                    {" << "\n";
        output << "                                        MDarray[ pos ] = read[ 'MDtag' ][j];" << "\n";
        output << "                                        pos = '';" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else pos += read[ 'MDtag' ][j];" << "\n";
        output << "                                }" << "\n";
        output << "                            });" << "\n";
        output << "" << "\n";
        output << "                            for( var j = 1; j < 8; ++j )" << "\n";
        output << "                                if( j in MDarray ) MDseed += ( j -1 ) + MDarray[j];" << "\n";
        output << "" << "\n";
        output << "                            var anno_split = seeds[i][ 'Name' ].split( '|' );" << "\n";
        output << "                            var annotation = anno_split[0];" << "\n";
        output << "" << "\n";
        output << "                            if( anno_split.length != 1 )" << "\n";
        output << "                            {" << "\n";
        output << "                                anno_split = anno_split[1].split( '_' );" << "\n";
        output << "                                annotation = annotation + '_' + anno_split[1]" << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            annotation = annotation" << "\n";
        output << "                                + ( MDseed == '' ? '' : ( '|' + MDseed ))" << "\n";
        output << "                                + ( RMSK == 'N' ? '' : '!' );" << "\n";
        output << "" << "\n";
        output << "                            $( '#seqs' ).append( \\\"<a id='seq\\\" + i + \\\"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\\\" + annotation + \\\"'></a>\\\" );" << "\n";
        output << "                            $( '#seq' + i + 'linked' ).css({" << "\n";
        output << "                                'text-decoration': 'none'," << "\n";
        output << "                                'color': 'black'" << "\n";
        output << "                                });" << "\n";
        output << "" << "\n";
        output << "                            $( '#seq' + i + 'linked' ).append( \\\"<div id='seq\\\" + i + \\\"' title='seqs\\\" + ( i + 1 ) + \\\"'>\\\" + json[ 'Sequence' ].charAt(i) + '</div>' );" << "\n";
        output << "                        }" << "\n";
        output << "                        else $( '#seqs' ).append( \\\"<div id='seq\\\" + i + \\\"' title='seqs\\\" + ( i + 1 ) + \\\"'>\\\" + json[ 'Sequence' ].charAt(i) + '</div>' );" << "\n";
        output << "" << "\n";
        output << "                        if( i != json[ 'Sequence' ].length-1 )" << "\n";
        output << "                            $( '#seq' + i ).css({ 'float': 'left' });" << "\n";
        output << "                    }\";" << "\n";
        output << "" << "\n";
        output << "                if( $isRNAfold ) echo \"" << "\n";
        output << "                    $( '#sequence' ).append( \\\"<div id='folds'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                    $( '#folds' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '30px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 0; i < json[ 'Fold0' ].length; ++i )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#folds' ).append( \\\"<div id='fold0\\\" + i + \\\"' title='seqs\\\" + ( i + 1 ) + \\\"'>\\\" + json[ 'Fold0' ].charAt(i) + '</div>' );" << "\n";
        output << "                        if( i != json[ 'Fold0' ].length-1 )" << "\n";
        output << "                            $( '#fold0' + i ).css({ 'float': 'left' });" << "\n";
        output << "                    }\";" << "\n";
        output << "" << "\n";
        output << "                echo \"" << "\n";
        output << "                    $( '#sequence' ).append( \\\"<div id='basecount2'></div>\\\" );" << "\n";
        output << "                    $( '#basecount2' ).css({ 'height': '10px' });" << "\n";
        output << "" << "\n";
        output << "                    for( var i = 1; i < json[ 'Sequence' ].length + extra; ++i ) if( i % 5 == 0 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#basecount2' ).append( \\\"<div id='base3-\\\" + i + \\\"' align='right'>\\\" + ( i == 5 ? '05' : i ) + '</div>' );" << "\n";
        output << "                        $( '#base3-' + i ).css({" << "\n";
        output << "                            'width': ( 5 * spc_num ) + 'px'," << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '10px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='seqalign'></div>\\\" );" << "\n";
        output << "                    $( '#seqalign' ).css({" << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#seqalign' ).append( \\\"<div id='reads'></div>\\\" );" << "\n";
        output << "                    $( '#reads' ).css({" << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra )* spc_num ) + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#seqalign' ).append( \\\"<div id='readslabel'></div>\\\" );" << "\n";
        output << "                    $( '#readslabel' ).css({" << "\n";
        output << "                        'width': '5px'," << "\n";
        output << "                        'height': ( datas.filter( Boolean ).length * read_height ) + 'px'," << "\n";
        output << "                        'background': 'linear-gradient( to top, $Color_Low 0%, $Color_Hight 100% )'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#readslabel' ).append( \\\"<div id='readslabeltop'>-\\\" + heat_max.toFixed(2) + '</div>' );" << "\n";
        output << "                    $( '#readslabeltop' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': seq_height + exp_height + shift_top1 + shift_top2 + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    if( datas.filter( Boolean ).length > 4 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#readslabel' ).append( \\\"<div id='readslabelq2'>-\\\" + ( heat_max * 0.75 ).toFixed(2) + '</div>' );" << "\n";
        output << "                        $( '#readslabelq2' ).css({" << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '15px'," << "\n";
        output << "                            'position': 'absolute'," << "\n";
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.25 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( datas.filter( Boolean ).length > 1 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#readslabel' ).append( \\\"<div id='readslabelmid'>-\\\" + ( heat_max * 0.5 ).toFixed(2) + '</div>' );" << "\n";
        output << "                        $( '#readslabelmid' ).css({" << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '15px'," << "\n";
        output << "                            'position': 'absolute'," << "\n";
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.5 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    if( datas.filter( Boolean ).length > 4 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#readslabel' ).append( \\\"<div id='readslabelq1'>-\\\" + ( heat_max * 0.25 ).toFixed(2) + '</div>' );" << "\n";
        output << "                        $( '#readslabelq1' ).css({" << "\n";
        output << "                            'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                            'font-size': '15px'," << "\n";
        output << "                            'position': 'absolute'," << "\n";
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.75 ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $( '#readslabel' ).append( \\\"<div id='readslabeldown'>-0.00</div>\\\" );" << "\n";
        output << "                    $( '#readslabeldown' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': ( datas.filter( Boolean ).length * read_height ) + seq_height + exp_height + shift_top1 + shift_top2 + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });\";" << "\n";
        output << "            }" << "\n";
        output << "            else echo \"" << "\n";
        output << "                    var reads = svg.append( 'g' )" << "\n";
        output << "                        .attr( 'id', 'reads' )" << "\n";
        output << "                        .attr( 'transform', 'translate( 26,' + ( fold_height + 180 ) + ')' )" << "\n";
        output << "                        .attr( 'font-family', 'Arial' )" << "\n";
        output << "                        .attr( 'font-size', '12px' )" << "\n";
        output << "                        ;" << "\n";
        output << "            \";" << "\n";
        output << "" << "\n";
        output << "            echo \"" << "\n";
        output << "                    var last_idx = -1;" << "\n";
        output << "                    var reads_count = 0;" << "\n";
        output << "" << "\n";
        output << "                    $.each( json[ 'Reads' ], function( idx, read )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( read[ 'Filter' ] == 'Y' ) return;" << "\n";
        output << "                        reads_count++;" << "\n";
        output << "" << "\n";
        output << "                        var pos = '';" << "\n";
        output << "                        var MDseed = '';" << "\n";
        output << "                        var MDarray = Array();" << "\n";
        output << "                        var TCarray = Array();" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'MDtag' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( !$.isNumeric( read[ 'MDtag' ][i] ))" << "\n";
        output << "                            {" << "\n";
        output << "                                MDarray[ pos ] = read[ 'MDtag' ][i];" << "\n";
        output << "                                pos = '';" << "\n";
        output << "                            }" << "\n";
        output << "                            else pos += read[ 'MDtag' ][i];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'TCtag' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( !$.isNumeric( read[ 'TCtag' ][i] ))" << "\n";
        output << "                            {" << "\n";
        output << "                                TCarray[ pos ] = read[ 'TCtag' ][i];" << "\n";
        output << "                                pos = '';" << "\n";
        output << "                            }" << "\n";
        output << "                            else pos += read[ 'TCtag' ][i];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 1; i < 8; ++i )" << "\n";
        output << "                            if( i in MDarray ) MDseed += ( i -1 ) + MDarray[i];\";" << "\n";
        output << "" << "\n";
        output << "            if( !$isSVG ) echo \"" << "\n";
        output << "                        if( last_idx != read[ 'Index' ] )" << "\n";
        output << "                        {" << "\n";
        output << "                            last_idx = read[ 'Index' ];" << "\n";
        output << "                            $( '#reads' ).append( \\\"<div id='selectedseed\\\" + last_idx + \\\"'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                            $( '#selectedseed' + last_idx ).css({" << "\n";
        output << "                                'position': 'relative'," << "\n";
        output << "                                'left': last_idx * spc_num + 'px'," << "\n";
        output << "                                'width': seeds[ last_idx + 4 ][ 'Length' ] * spc_num + 'px'," << "\n";
        output << "                                });" << "\n";
        output << "" << "\n";
        output << "                            var anno_split = seeds[ last_idx + 4 ][ 'Name' ].split( '|' );" << "\n";
        output << "                            var annotation = anno_split[0];" << "\n";
        output << "" << "\n";
        output << "                            if( anno_split.length != 1 )" << "\n";
        output << "                            {" << "\n";
        output << "                                anno_split = anno_split[1].split( '_' );" << "\n";
        output << "                                annotation = annotation + '_' + anno_split[1]" << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            $( '#selectedseed' + last_idx ).append( \\\"<a id='read\\\" + last_idx + \\\"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\\\"" << "\n";
        output << "                                + annotation" << "\n";
        output << "                                + ( MDseed == '' ? '' : ( '|' + MDseed ))" << "\n";
        output << "                                + ( read[ 'isRMSK' ] == 'N' ? '' : '!' )" << "\n";
        output << "                                + \\\"'></a>\\\" );" << "\n";
        output << "" << "\n";
        output << "                            $( '#read' + last_idx + 'linked' ).css({ 'text-decoration': 'none' });" << "\n";
        output << "                        }\";" << "\n";
        output << "" << "\n";
        output << "            if( $isSVG ) echo\"" << "\n";
        output << "                        for( var i = 0; i < read[ 'Length' ]; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( i in MDarray )" << "\n";
        output << "                            {" << "\n";
        output << "                                var md = reads.append( 'text' )" << "\n";
        output << "                                    .attr( 'x', ( read[ 'Index' ] + i ) * spc_num )" << "\n";
        output << "                                    .attr( 'y', reads_count * read_height )" << "\n";
        output << "                                    .attr( 'fill-opacity', '0.4' )" << "\n";
        output << "                                    .text( MDarray[i] )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                                if( MDarray[i] == 'A' ) md.attr( 'fill', colorA );" << "\n";
        output << "                                if( MDarray[i] == 'C' ) md.attr( 'fill', colorC );" << "\n";
        output << "                                if( MDarray[i] == 'G' ) md.attr( 'fill', colorG );" << "\n";
        output << "                                if( MDarray[i] == 'U' ) md.attr( 'fill', colorU );" << "\n";
        output << "                            }" << "\n";
        output << "                            else if( i in TCarray )" << "\n";
        output << "                            {" << "\n";
        output << "                                var tc = reads.append( 'text' )" << "\n";
        output << "                                    .attr( 'x', ( read[ 'Index' ] + i ) * spc_num )" << "\n";
        output << "                                    .attr( 'y', reads_count * read_height )" << "\n";
        output << "                                    .attr( 'fill-opacity', '0.4' )" << "\n";
        output << "                                    .text( TCarray[i] )" << "\n";
        output << "                                    ;" << "\n";
        output << "" << "\n";
        output << "                                if( TCarray[i] == 'A' ) tc.attr( 'fill', colorA );" << "\n";
        output << "                                if( TCarray[i] == 'C' ) tc.attr( 'fill', colorC );" << "\n";
        output << "                                if( TCarray[i] == 'G' ) tc.attr( 'fill', colorG );" << "\n";
        output << "                                if( TCarray[i] == 'U' ) tc.attr( 'fill', colorU );" << "\n";
        output << "                            }" << "\n";
        output << "                            else" << "\n";
        output << "                            {" << "\n";
        output << "                                reads.append( 'text' )" << "\n";
        output << "                                    .attr( 'x', ( read[ 'Index' ] + i ) * spc_num )" << "\n";
        output << "                                    .attr( 'y', reads_count * read_height )" << "\n";
        output << "                                    .attr( 'fill', color_map( heat_array[ idx ]) )" << "\n";
        output << "                                    .text( json[ 'Sequence' ][ read[ 'Index' ] + i ] )" << "\n";
        output << "                                    ;" << "\n";
        output << "                            }" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'Tail' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            var tail = reads.append( 'text' )" << "\n";
        output << "                                .attr( 'x', ( read[ 'Index' ] + read[ 'Length' ] + i ) * spc_num )" << "\n";
        output << "                                .attr( 'y', reads_count * read_height )" << "\n";
        output << "                                .text( read[ 'Tail' ].charAt(i) )" << "\n";
        output << "                                ;" << "\n";
        output << "" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'A' ) tail.attr( 'fill', colorA );" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'C' ) tail.attr( 'fill', colorC );" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'G' ) tail.attr( 'fill', colorG );" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'U' ) tail.attr( 'fill', colorU );" << "\n";
        output << "                        }\";" << "\n";
        output << "" << "\n";
        output << "            else echo\"" << "\n";
        output << "                        $( '#read' + last_idx + 'linked' ).append( \\\"<div id='read\\\" + idx + \\\"' title=Read\\\" + idx + \\\"></div>\\\" );" << "\n";
        output << "                        $( '#read' + idx ).css({" << "\n";
        output << "                            'width': (( read[ 'Length' ] + read[ 'Tail' ].length ) * spc_num ) + 'px'," << "\n";
        output << "                            });" << "\n";
        output << "" << "\n";
        output << "                        pos = 0;" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'Length' ]; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( i in MDarray )" << "\n";
        output << "                            {" << "\n";
        output << "                                $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'gm' + i + \\\"'></div>\\\" );" << "\n";
        output << "                                $( '#read' + idx + 'gm' + i ).css({" << "\n";
        output << "                                    'width': ( pos * spc_num ) + 'px'," << "\n";
        output << "                                    'height': read_height + 'px'," << "\n";
        output << "                                    'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                    'position': 'relative'," << "\n";
        output << "                                    'float': 'left'," << "\n";
        output << "                                    'z-index': '-1'" << "\n";
        output << "                                    });" << "\n";
        output << "" << "\n";
        output << "                                pos = 0;" << "\n";
        output << "" << "\n";
        output << "                                $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'md' + i + \\\"'>\\\" + MDarray[i] + '</div>' );" << "\n";
        output << "                                $( '#read' + idx + 'md' + i ).css({" << "\n";
        output << "                                    'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                                    'font-weight': '100'," << "\n";
        output << "                                    'font-size': '18px'," << "\n";
        output << "                                    'text-align': 'center'," << "\n";
        output << "                                    'height': read_height + 'px'," << "\n";
        output << "                                    'width': spc_num + 'px'," << "\n";
        output << "                                    'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                    'color': 'lightpink'," << "\n";
        output << "                                    'position': 'relative'," << "\n";
        output << "                                    'float': 'left'," << "\n";
        output << "                                    'z-index': '-1'" << "\n";
        output << "                                    });" << "\n";
        output << "                            }" << "\n";
        output << "                            else if( i in TCarray )" << "\n";
        output << "                            {" << "\n";
        output << "                                $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'gm' + i + \\\"'></div>\\\" );" << "\n";
        output << "                                $( '#read' + idx + 'gm' + i ).css({" << "\n";
        output << "                                    'width': ( pos * spc_num ) + 'px'," << "\n";
        output << "                                    'height': read_height + 'px'," << "\n";
        output << "                                    'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                    'position': 'relative'," << "\n";
        output << "                                    'float': 'left'," << "\n";
        output << "                                    'z-index': '-1'" << "\n";
        output << "                                    });" << "\n";
        output << "" << "\n";
        output << "                                pos = 0;" << "\n";
        output << "" << "\n";
        output << "                                $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'tc' + i + \\\"'>\\\" + TCarray[i] + '</div>' );" << "\n";
        output << "                                $( '#read' + idx + 'tc' + i ).css({" << "\n";
        output << "                                    'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                                    'font-weight': '100'," << "\n";
        output << "                                    'font-size': '18px'," << "\n";
        output << "                                    'text-align': 'center'," << "\n";
        output << "                                    'height': read_height + 'px'," << "\n";
        output << "                                    'width': spc_num + 'px'," << "\n";
        output << "                                    'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                    'color': 'lightblue'," << "\n";
        output << "                                    'position': 'relative'," << "\n";
        output << "                                    'float': 'left'," << "\n";
        output << "                                    'z-index': '-1'" << "\n";
        output << "                                    });" << "\n";
        output << "                            }" << "\n";
        output << "                            else" << "\n";
        output << "                            {" << "\n";
        output << "                                pos++;" << "\n";
        output << "" << "\n";
        output << "                                if( i == read[ 'Length' ] -1 )" << "\n";
        output << "                                {" << "\n";
        output << "                                    $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'gm' + i + \\\"'></div>\\\" );" << "\n";
        output << "                                    $( '#read' + idx + 'gm' + i ).css({" << "\n";
        output << "                                        'width': ( pos * spc_num ) + 'px'," << "\n";
        output << "                                        'height': read_height + 'px'," << "\n";
        output << "                                        'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                        'position': 'relative'," << "\n";
        output << "                                        'float': 'left'," << "\n";
        output << "                                        'z-index': '-1'" << "\n";
        output << "                                        });" << "\n";
        output << "" << "\n";
        output << "                                    if( read[ 'Tail' ] == '' )" << "\n";
        output << "                                    {" << "\n";
        output << "                                        $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'gmEnd' + \\\"'></div>\\\" );" << "\n";
        output << "                                        $( '#read' + idx + 'gmEnd' ).css({" << "\n";
        output << "                                            'width': ( 0 * spc_num ) + 'px'," << "\n";
        output << "                                            'height': read_height + 'px'," << "\n";
        output << "                                            'background': color_map( heat_array[ idx ])," << "\n";
        output << "                                            'position': 'relative'," << "\n";
        output << "                                            'z-index': '-1'" << "\n";
        output << "                                            });" << "\n";
        output << "                                    }" << "\n";
        output << "                                }" << "\n";
        output << "                            }" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'Tail' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            $( '#read' + idx ).append( \\\"<div id='read\\\" + idx + 'pm' + i + \\\"'>\\\" + read[ 'Tail' ].charAt(i) + '</div>' );" << "\n";
        output << "                            $( '#read' + idx + 'pm' + i ).css({" << "\n";
        output << "                                'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                                'font-size': '18px'," << "\n";
        output << "                                'text-align': 'center'," << "\n";
        output << "                                'height': read_height + 'px'," << "\n";
        output << "                                'position': 'relative'," << "\n";
        output << "                                'z-index': '-1'" << "\n";
        output << "                                });" << "\n";
        output << "" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'A' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorA });" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'C' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorC });" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'G' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorG });" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'U' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorU });" << "\n";
        output << "" << "\n";
        output << "                            if( i != read[ 'Tail' ].length -1 )" << "\n";
        output << "                            {" << "\n";
        output << "                                $( '#read' + idx + 'pm' + i ).css({" << "\n";
        output << "                                    'width': spc_num + 'px'," << "\n";
        output << "                                    'float': 'left'" << "\n";
        output << "                                    });" << "\n";
        output << "                            }" << "\n";
        output << "                        }\";" << "\n";
        output << "" << "\n";
        output << "            echo \"                    });\";" << "\n";
        output << "" << "\n";
        output << "            if( $isSVG ) echo\"" << "\n";
        output << "                    svg.attr( 'height', ( reads_count * read_height + fold_height + 180 ) + 'px' );" << "\n";
        output << "                \";" << "\n";
        output << "            else echo\"" << "\n";
        output << "                    $( '#chart' ).css({" << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        'height': reads_count * read_height + seq_height + exp_height + 'px'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='line'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                    $( '#line' ).css({" << "\n";
        output << "                        'margin-top': fold_height + 'px'," << "\n";
        output << "                        'height': $( '#chart' ).height() + shift_top2 + 'px'," << "\n";
        output << "                        'width': '1px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': '40px'," << "\n";
        output << "                        'left': '-1px'," << "\n";
        output << "                        'background-color': 'indianred'," << "\n";
        output << "                        'pointer-events': 'none'" << "\n";
        output << "                        });\";" << "\n";
        output << "" << "\n";
        output << "            echo \"                });\";" << "\n";
        output << "" << "\n";
        output << "            if( !$isSVG )" << "\n";
        output << "            {" << "\n";
        output << "                echo \"" << "\n";
        output << "                var selectedfold = [];" << "\n";
        output << "                var selectedread = $( '#chart' );" << "\n";
        output << "                var selectedseed = $( '#chart' );" << "\n";
        output << "" << "\n";
        output << "                $( document ).mousemove( function( event )" << "\n";
        output << "                {" << "\n";
        output << "                    if( event.pageY > shift_top1 )" << "\n";
        output << "                    {" << "\n";
        output << "                        $( '#line' ).css({" << "\n";
        output << "                            'left': ( event.pageX > $( '#chart' ).width() ? -1 : event.pageX ) + 'px'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "                $( function()" << "\n";
        output << "                {" << "\n";
        output << "                    $( document ).tooltip(" << "\n";
        output << "                    {" << "\n";
        output << "                        content: function()" << "\n";
        output << "                        {" << "\n";
        output << "                            var tip = '';" << "\n";
        output << "                            var data = Array();" << "\n";
        output << "                            selectedread = $( this );" << "\n";
        output << "" << "\n";
        output << "                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'fold' )" << "\n";
        output << "                                return tip;" << "\n";
        output << "" << "\n";
        output << "                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'expr' )" << "\n";
        output << "                            {" << "\n";
        output << "                                data = exprs[ selectedread.attr( 'title' ).substr( 4 ) ];" << "\n";
        output << "                                tip += 'PPM: ' + ( data[ 'PPM' ] * expr_max / 100 ).toFixed(2);\";" << "\n";
        output << "" << "\n";
        output << "                if( $isRNAfold ) echo \"" << "\n";
        output << "                                tip += '</br>Entropy: ' + ( etpys[ selectedread.attr( 'title' ).substr( 4 ) ] ).toFixed(2);\";" << "\n";
        output << "                echo \"" << "\n";
        output << "                                if( data[ 'A' ] != 0 ) tip += '</br>A: ' + ( data[ 'A' ] * 100 ).toFixed(3) + '%';" << "\n";
        output << "                                if( data[ 'C' ] != 0 ) tip += '</br>C: ' + ( data[ 'C' ] * 100 ).toFixed(3) + '%';" << "\n";
        output << "                                if( data[ 'G' ] != 0 ) tip += '</br>G: ' + ( data[ 'G' ] * 100 ).toFixed(3) + '%';" << "\n";
        output << "                                if( data[ 'U' ] != 0 ) tip += '</br>T: ' + ( data[ 'U' ] * 100 ).toFixed(3) + '%';" << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )" << "\n";
        output << "                            {" << "\n";
        output << "                                data = datas[ selectedread.attr( 'title' ).substr( 4 ) ];" << "\n";
        output << "" << "\n";
        output << "                                var pos = '';" << "\n";
        output << "                                var MDseed = '';" << "\n";
        output << "                                var seq = sequ.substr( data[ 'Index' ], data[ 'Length' ]);" << "\n";
        output << "                                var MDarray = Array();" << "\n";
        output << "                                var TCarray = Array();" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 0; i < data[ 'MDtag' ].length; ++i )" << "\n";
        output << "                                {" << "\n";
        output << "                                    if( !$.isNumeric( data[ 'MDtag' ][i] ))" << "\n";
        output << "                                    {" << "\n";
        output << "                                        MDarray[ pos ] = data[ 'MDtag' ][i];" << "\n";
        output << "                                        pos = '';" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else pos += data[ 'MDtag' ][i];" << "\n";
        output << "                                }" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 0; i < data[ 'TCtag' ].length; ++i )" << "\n";
        output << "                                {" << "\n";
        output << "                                    if( !$.isNumeric( data[ 'TCtag' ][i] ))" << "\n";
        output << "                                    {" << "\n";
        output << "                                        TCarray[ pos ] = data[ 'TCtag' ][i];" << "\n";
        output << "                                        pos = '';" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else pos += data[ 'TCtag' ][i];" << "\n";
        output << "                                }" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 1; i < 8; ++i )" << "\n";
        output << "                                    if( i in MDarray ) MDseed += ( i -1 ) + MDarray[i];" << "\n";
        output << "" << "\n";
        output << "                                $( '#selectbox' ).css({" << "\n";
        output << "                                    'left': ( data[ 'Index' ] * spc_num + 4 ) + 'px'," << "\n";
        output << "                                    'width': (( data[ 'Length' ] ) * spc_num ) + 'px'," << "\n";
        output << "                                    'display': 'block'" << "\n";
        output << "                                    });" << "\n";
        output << "" << "\n";
        output << "                                selectedfold = [];" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 0; i < data[ 'Length' ]; ++i )" << "\n";
        output << "                                {" << "\n";
        output << "                                    selectedfold[i] = '#fold' + ( data[ 'Index' ] +i +1 );" << "\n";
        output << "                                    $( selectedfold[i] ).css({" << "\n";
        output << "                                        'background-color': 'Gold'" << "\n";
        output << "                                        });" << "\n";
        output << "                                }" << "\n";
        output << "" << "\n";
        output << "                                tip += 'Annotation: ' + anno;" << "\n";
        output << "                                tip += \".( $Data_Array[0][5] != '.'? ( \"'-' + data[ 'Arm' ]\" ) : \"''\" ).\";" << "\n";
        output << "" << "\n";
        output << "                                tip += '_' + sequ.substr(( data[ 'Index' ] + 1 ), 7 );" << "\n";
        output << "                                tip += ( MDseed == '' ? '' : ( '|<b style=\\\"color:red;\\\">' + MDseed + '</b>' ));" << "\n";
        output << "" << "\n";
        output << "                                tip += ( data[ 'isRMSK' ] == 'N' ? '' : ' (RMSK)' ) + '</br>';" << "\n";
        output << "                                tip += 'Sequence: ';" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 0; i < data[ 'Length' ]; ++i )" << "\n";
        output << "                                {" << "\n";
        output << "                                    if( i in MDarray )" << "\n";
        output << "                                    {" << "\n";
        output << "                                        tip += '<b style=\\\"color:red;\\\">' + MDarray[i] + '</b>';" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else if( i in TCarray )" << "\n";
        output << "                                    {" << "\n";
        output << "                                        tip += '<b style=\\\"color:blue;\\\">' + TCarray[i] + '</b>';" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else tip += seq[i];" << "\n";
        output << "                                }" << "\n";
        output << "" << "\n";
        output << "                                tip += '</br>';" << "\n";
        output << "" << "\n";
        output << "                                if( data[ 'MDtag' ] != '' ) tip += '<font color=red>MDtag: ' + data[ 'MDtag' ] + '</font></br>'; " << "\n";
        output << "                                if( data[ 'TCtag' ] != '' ) tip += '<font color=blue>TCtag: ' + data[ 'TCtag' ] + '</font></br>'; " << "\n";
        output << "" << "\n";
        output << "                                tip += 'Length: ' + data[ 'Length' ];" << "\n";
        output << "                                tip += ( data[ 'Tail' ] == '' ? '</br>'" << "\n";
        output << "                                     : ( ' + ' + data[ 'Tail' ].length + '</br>Tail: ' + data[ 'Tail' ] + '</br>' ));" << "\n";
        output << "" << "\n";
        output << "                                tip += 'PPM: ' + data[ 'PPM' ].toFixed(2); " << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'seqs' )" << "\n";
        output << "                            {" << "\n";
        output << "                                selectedseed = $( '#selectedseed' + ( selectedread.attr( 'title' ).substr( 4 ) - 5 ));" << "\n";
        output << "" << "\n";
        output << "                                var idx = (( selectedread.attr( 'title' ).substr( 4 ) < 5 ) ? 0 :" << "\n";
        output << "                                          (( selectedread.attr( 'title' ).substr( 4 ) > sequ.length -3 ) ? ( sequ.length -7 ) :" << "\n";
        output << "                                           ( selectedread.attr( 'title' ).substr( 4 ) - 4 )));" << "\n";
        output << "" << "\n";
        output << "                                $( '#selectbox' ).css({" << "\n";
        output << "                                    'left': ( idx + 0.225 ) * spc_num + 'px'," << "\n";
        output << "                                    'width': ( 7 * spc_num ) + 'px'," << "\n";
        output << "                                    'display': 'block'" << "\n";
        output << "                                    });" << "\n";
        output << "" << "\n";
        output << "                                selectedseed.css({" << "\n";
        output << "                                    'border': '2px solid gold'," << "\n";
        output << "                                    'border-radius': '12px'," << "\n";
        output << "                                    });" << "\n";
        output << "" << "\n";
        output << "                                if(( selectedread.attr( 'title' ).substr( 4 ) - 1 ) in seeds )" << "\n";
        output << "                                {" << "\n";
        output << "                                    tip += seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][ 'Name' ];" << "\n";
        output << "                                    tip += '</br>GMPM: ' + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][ 'GMPM' ].toFixed(2);" << "\n";
        output << "                                    tip += '</br>GM: '   + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][  'GM'  ].toFixed(2);" << "\n";
        output << "                                    tip += '</br>PM: '   + seeds[ selectedread.attr( 'title' ).substr( 4 ) - 1 ][  'PM'  ].toFixed(2);" << "\n";
        output << "                                }" << "\n";
        output << "" << "\n";
        output << "                                return tip;" << "\n";
        output << "                            }" << "\n";
        output << "" << "\n";
        output << "                            selectedread.css({" << "\n";
        output << "                                'border': '2px solid gold'," << "\n";
        output << "                                'border-radius': '12px'" << "\n";
        output << "                                });" << "\n";
        output << "" << "\n";
        output << "                            return tip;" << "\n";
        output << "                        }," << "\n";
        output << "" << "\n";
        output << "                        close: function( event, ui )" << "\n";
        output << "                        {" << "\n";
        output << "                            ui.tooltip.hover(" << "\n";
        output << "                                function ()" << "\n";
        output << "                                {" << "\n";
        output << "                                    $( this ).stop( true ).fadeTo( 400, 1 ); " << "\n";
        output << "" << "\n";
        output << "                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )" << "\n";
        output << "                                        $( '#selectbox' ).css({ 'display': 'block' });" << "\n";
        output << "" << "\n";
        output << "                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'seqs' )" << "\n";
        output << "                                    {" << "\n";
        output << "                                        var idx = (( selectedread.attr( 'title' ).substr( 4 ) < 5 ) ? 0 :" << "\n";
        output << "                                                  (( selectedread.attr( 'title' ).substr( 4 ) > sequ.length -3 ) ? ( sequ.length -7 ) :" << "\n";
        output << "                                                   ( selectedread.attr( 'title' ).substr( 4 ) - 4 )));" << "\n";
        output << "" << "\n";
        output << "                                        $( '#selectbox' ).css({" << "\n";
        output << "                                            'left': idx * spc_num + 'px'," << "\n";
        output << "                                            'width': 7 * spc_num + 'px'," << "\n";
        output << "                                            'display': 'block'" << "\n";
        output << "                                            });" << "\n";
        output << "" << "\n";
        output << "                                        selectedseed.css({" << "\n";
        output << "                                            'border': '2px solid gold'," << "\n";
        output << "                                            'border-radius': '12px'" << "\n";
        output << "                                            });" << "\n";
        output << "                                    }" << "\n";
        output << "                                    else" << "\n";
        output << "                                    {" << "\n";
        output << "                                        selectedread.css({" << "\n";
        output << "                                            'border': '2px solid gold'," << "\n";
        output << "                                            'border-radius': '12px'" << "\n";
        output << "                                            });" << "\n";
        output << "                                    }" << "\n";
        output << "                                }," << "\n";
        output << "                                function ()" << "\n";
        output << "                                {" << "\n";
        output << "                                    $( this ).fadeOut( '400', function(){ $( this ).remove(); });" << "\n";
        output << "" << "\n";
        output << "                                    if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'Read' )" << "\n";
        output << "                                        $( '#selectbox' ).css({ 'display': 'none' });" << "\n";
        output << "" << "\n";
        output << "                                    selectedread.css({ 'border': 'none' });" << "\n";
        output << "                                    selectedseed.css({ 'border': 'none' });" << "\n";
        output << "                                }" << "\n";
        output << "                            );" << "\n";
        output << "" << "\n";
        output << "                            $( '#selectbox' ).css({ 'display': 'none' });" << "\n";
        output << "                            selectedread.css({ 'border': 'none' });" << "\n";
        output << "                            selectedseed.css({ 'border': 'none' });" << "\n";
        output << "" << "\n";
        output << "                            for( var i = 0; i < selectedfold.length; ++i )" << "\n";
        output << "                            {" << "\n";
        output << "                                $( selectedfold[i] ).css({" << "\n";
        output << "                                    'background-color': 'transparent'" << "\n";
        output << "                                    });" << "\n";
        output << "                            }" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "                });\";" << "\n";
        output << "            }" << "\n";
        output << "" << "\n";
        output << "            echo\"" << "\n";
        output << "                </script>\";" << "\n";
        output << "        }" << "\n";
        output << "    ?>" << "\n";
        output << "    </body>" << "\n";
        output << "</html>" << "\n";

        output.close();
    }
};

} // end of namespace algorithm
} // end of namespace ago
