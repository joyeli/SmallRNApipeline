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
            rawbed.getReadSeq( *genome ).substr( 1, 7 ),
            get_arm( rawbed ),
            ( rawbed.getTail() != "" ? rawbed.getTail() : "." )
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
                    case 'A': c = 'T'; break;
                    case 'T': c = 'A'; break;
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
            std::size_t& extend_refseq
            )
    {
        std::ofstream outidx, outtsv;
        std::map< std::string, SeqType > chr_mapping;
        std::size_t idx;

        for( std::size_t smp = 0; smp < bed_samples.size(); ++smp )
        {
            idx = 0;
            chr_mapping = get_chrmap_table( bed_samples[ smp ].second, biotype, genome_table, extend_refseq );

            outidx.open( output_name + bed_samples[ smp ].first + ".idx" );
            outtsv.open( output_name + bed_samples[ smp ].first + ".tsv" );

            for( auto& mir : chr_mapping )
            {
                outidx << mir.first << "\t" << idx << "\n";;
                outtsv << mir.first << mir.second;
                idx += mir.second.line_counts();
            }

            outidx.close();
            outtsv.close();
        }
    }

    static std::map< std::string, SeqType > get_chrmap_table(
            // std::vector< AnnotationRawBed<> >& smp_anno,
            std::vector< ago::format::MDRawBed >& smp_anno,
            const std::string& biotype,
            auto& genome_table,
            std::size_t& extend_refseq
            )
    {
        std::map< std::string, SeqType > chr_mappings;
        std::string anno;

        for( auto& raw_bed : smp_anno )
        {
            // for( std::size_t i = 0; i < raw_bed.annotation_info_.size(); ++i )
            {
                std::size_t i = 0; // do first priority
                if( i < raw_bed.annotation_info_.size() && !raw_bed.annotation_info_[i].empty() )
                {
                    if(  raw_bed.start_ < extend_refseq ) continue;
                    if(( raw_bed.annotation_info_[i][0] == biotype ) ||
                       ( biotype == "miRNA_mirtron" && ( raw_bed.annotation_info_[i][0] == "miRNA" || raw_bed.annotation_info_[i][0] == "mirtron" ))) 
                    {
                        for( std::size_t j = 0; j < raw_bed.annotation_info_[i].size(); j+=2 )
                        {
                            anno = biotype.substr( 0, 5 ) != "miRNA" && biotype != "mirtron"
                                 ? raw_bed.annotation_info_[i][ j+1 ]
                                 : raw_bed.annotation_info_[i][ j+1 ].substr( 0, raw_bed.annotation_info_[i][ j+1 ].length() -3 );

                            if( chr_mappings.find( anno ) == chr_mappings.end() )
                            {
                                chr_mappings[ anno ] = SeqType();
                                chr_mappings[ anno ].init( raw_bed, genome_table, raw_bed.annotation_info_[i][0] );
                            }

                            chr_mappings[ anno ].insert( raw_bed );
                        }
                    }
                }
            }
        }

        for( auto& mir : chr_mappings ) mir.second.formation( extend_refseq );
        return chr_mappings;
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
        output << "" << "\n";
        output << "        if( $_GET['TSV_File'] == '' ) $TSV_File = $_POST['TSV_File'];" << "\n";
        output << "        if( $_GET['Annotation_Select'] == '' ) $Annotation_Select = $_POST['Annotation_Select'];" << "\n";
        output << "" << "\n";
        output << "        $Length = $_POST['Length'];" << "\n";
        output << "        $Filter = $_POST['Filter'];" << "\n";
        output << "        $Color_Hight = $_POST['Color_Hight'];" << "\n";
        output << "        $Color_Low = $_POST['Color_Low'];" << "\n";
        output << "        $Segment_Select = $_POST['Segment_Select'];" << "\n";
        output << "" << "\n";
        output << "        if( $Color_Hight == '' ) $Color_Hight = 'Black';" << "\n";
        output << "        if( $Color_Low == '' ) $Color_Low = 'WhiteSmoke';" << "\n";
        output << "" << "\n";
        output << "#<!--================== TSV File ====================-->" << "\n";
        output << "                " << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "" << "\n";
        output << "        $TSV = Shell_Exec( 'ls | grep .tsv' );" << "\n";
        output << "        $TSV_List = Explode( \"\\n\", $TSV );" << "\n";
        output << "        $List_Size = Count( $TSV_List );" << "\n";
        output << "" << "\n";
        output << "        echo '<select name=TSV_File onchange=this.form.submit();>';" << "\n";
        output << "        echo '<option '; if($TSV_File=='') echo 'selected'; echo '>Select TSV</option>';" << "\n";
        output << "" << "\n";
        output << "        For( $i = 0; $i < $List_Size-1; ++$i )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$TSV_List[$i].' ';" << "\n";
        output << "" << "\n";
        output << "            if( $TSV_File == $TSV_List[$i] ) " << "\n";
        output << "                echo 'selected ';" << "\n";
        output << "" << "\n";
        output << "            echo '>'.$TSV_List[$i].'</option>';" << "\n";
        output << "        }" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Segment_Select' value='$Segment_Select' />" << "\n";
        output << "            </form>\";" << "\n";
        output << "" << "\n";
        output << "#<!--================== Annotation Select ====================-->" << "\n";
        output << "" << "\n";
        output << "        $inFile = File_get_contents( substr( $TSV_File, 0, strlen( $TSV_File ) -4 ).'.idx' );" << "\n";
        output << "        $inFile_Lines = Explode( \"\\n\", $inFile );" << "\n";
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
        output << "            if( $Annotation_Select == $anno[0] && $Annotation_Select != '' )" << "\n";
        output << "            {" << "\n";
        output << "                $TSV = new SplFileObject( $TSV_File );" << "\n";
        output << "                $TSV->seek( $anno[1] );" << "\n";
        output << "" << "\n";
        output << "                $annoSeq = Explode( \"\\t\", trim( $TSV->current(), \"\\n\" ));" << "\n";
        output << "                $Full_Sequc = $annoSeq[1];" << "\n";
        output << "                $Ref_Chr    = $annoSeq[2];" << "\n";
        output << "                $Ref_Start  = $annoSeq[3];" << "\n";
        output << "                $Ref_Strand = $annoSeq[4];" << "\n";
        output << "" << "\n";
        output << "                while( true )" << "\n";
        output << "                {" << "\n";
        output << "                    $TSV->next();" << "\n";
        output << "                    $tsvData = Explode( \"\\t\", trim( $TSV->current(), \"\\n\" ));" << "\n";
        output << "                    if( $tsvData[0] != '' ) break;" << "\n";
        output << "                    Array_Push( $Data_Array, $tsvData );" << "\n";
        output << "                }" << "\n";
        output << "" << "\n";
        output << "                $TSV = null;" << "\n";
        output << "            }" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        Ksort( $Anno_Array );" << "\n";
        output << "" << "\n";
        output << "        echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "        echo '<select name=Annotation_Select onchange=this.form.submit();>';" << "\n";
        output << "        echo \"<option value='' \"; if($Annotation_Select=='') echo 'selected'; echo '>Select Annotations</option>';" << "\n";
        output << "" << "\n";
        output << "        Foreach( $Anno_Array as $anno => $idx )" << "\n";
        output << "        {" << "\n";
        output << "            echo '<option value='.$anno.' ';" << "\n";
        output << "            if( $Annotation_Select == $anno )  echo 'selected ';" << "\n";
        output << "            echo '>'.$anno.'</option>';" << "\n";
        output << "        }" << "\n";
        output << "" << "\n";
        output << "        echo \"</select>" << "\n";
        output << "            <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "            <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "            <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "            <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            </form>\";" << "\n";
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
        output << "            for( $j = Count( $Segm_Array ); $j < Count( $Data_Array ); ++$j )" << "\n";
        output << "                $Segm_Array[$j] = $Segm_Range;" << "\n";
        output << "" << "\n";
        output << "            echo '<form action='.$_SERVER['PHP_SELF'].' method=post style=display:inline;>';" << "\n";
        output << "            echo '<select name=Segment_Select onchange=this.form.submit();>';" << "\n";
        output << "            echo \"<option value='' \"; if($Segment_Select=='') echo 'selected'; echo '>Select Segments</option>';" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )" << "\n";
        output << "            {" << "\n";
        output << "                echo '<option value=\"'.$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).'ppm)\" ';" << "\n";
        output << "" << "\n";
        output << "                if( $Segment_Select == $Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).'ppm)' )" << "\n";
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
        output << "                <input type='hidden' name='Length' value='$Length' />" << "\n";
        output << "                <input type='hidden' name='Filter' value='$Filter' />" << "\n";
        output << "                <input type='hidden' name='Color_Hight' value='$Color_Hight' />" << "\n";
        output << "                <input type='hidden' name='Color_Low' value='$Color_Low' />" << "\n";
        output << "                <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "                <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "                </form>\";" << "\n";
        output << "" << "\n";
        output << "            echo \"<script>var select_color_map = d3.scaleLinear().domain([ 0, 100 ]).range([ '$Color_Low', '$Color_Hight' ]);\";" << "\n";
        output << "" << "\n";
        output << "            for( $i = 0; $i < Count( $Segm_Uniqs ); ++$i )" << "\n";
        output << "                echo \"$( 'option[value=\\\"\".$Segm_Uniqs[$i][0].'-'.$Segm_Uniqs[$i][1].' ('.Number_Format( (float)$Segm_Uniqs[$i][2], 2, '.', '' ).\"ppm)\\\"]' ).css(" << "\n";
        output << "                {" << "\n";
        output << "                    'background-color': select_color_map( '\".( $Segm_Uniqs[$i][2] * $Total_PPM ).\"' )," << "\n";
        output << "                    'color': '\".(( $Segm_Uniqs[$i][2] * $Total_PPM ) > 25 ? 'white' : 'black'  ).\"'" << "\n";
        output << "                });\";" << "\n";
        output << "" << "\n";
        output << "            echo '</script>';" << "\n";
        output << "        }" << "\n";
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
        output << "#<!--==================== Heatmap Color ======================-->" << "\n";
        output << "" << "\n";
        output << "        echo '<input type=text name=Color_Low size=8 value='.$Color_Low.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "        echo '<input type=text name=Color_Hight size=8 value='.$Color_Hight.' onfocus=\"{this.value=\\'\\';}\" />';" << "\n";
        output << "" << "\n";
        output << "        echo \"" << "\n";
        output << "            <input type='hidden' name='TSV_File' value='$TSV_File' />" << "\n";
        output << "            <input type='hidden' name='Annotation_Select' value='$Annotation_Select' />" << "\n";
        output << "            <input type='hidden' name='Segment_Select' value='$Segment_Select' />" << "\n";
        output << "            <input type='submit' value='Submit' /> " << "\n";
        output << "            </form>\";" << "\n";
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
        output << "                    Fwrite( $Ftemp, \"      \\\"Filter\\\" : \\\"Y\\\"\\n\" );" << "\n";
        output << "                else" << "\n";
        output << "                    Fwrite( $Ftemp, \"      \\\"Filter\\\" : \\\"N\\\"\\n\" );" << "\n";
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
        output << "                var seeds = Array();" << "\n";
        output << "" << "\n";
        output << "                var anno = '';" << "\n";
        output << "                var sequ = '';" << "\n";
        output << "" << "\n";
        output << "                var colorA = 'red';" << "\n";
        output << "                var colorC = 'blue';" << "\n";
        output << "                var colorG = 'goldenrod';" << "\n";
        output << "                var colorT = 'green';" << "\n";
        output << "                var colorP = 'darkgray';" << "\n";
        output << "" << "\n";
        output << "                var exp_height = 120;" << "\n";
        output << "                var seq_height = 50;" << "\n";
        output << "                var read_height = 20;" << "\n";
        output << "                var spc_num = 18;" << "\n";
        output << "                var expr_max = 0;" << "\n";
        output << "                var shift_top = 33;" << "\n";
        output << "" << "\n";
        output << "                var seg_start = $Seg_Start;" << "\n";
        output << "                var seg_end   = $Seg_End;" << "\n";
        output << "" << "\n";
        output << "                $.getJSON( '$Temp', function( json )" << "\n";
        output << "                {" << "\n";
        output << "                    anno = json[ 'Annotation' ];" << "\n";
        output << "                    sequ = json[ 'Sequence'   ];" << "\n";
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
        output << "                        expr_array[ i ][  'T'  ] = 0;" << "\n";
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
        output << "" << "\n";
        output << "                        for( var i = start; i < ( end + read[ 'Tail' ].length ); ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            expr_array[ i ][ 'PPM' ] += read[ 'PPM' ];" << "\n";
        output << "" << "\n";
        output << "                            if( i < end ) expr_array[ i ][  'P'  ] += read[ 'PPM' ];" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
        output << "                        for( var i = 0; i < read[ 'Tail' ].length; ++i )" << "\n";
        output << "                        {" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'A' ) expr_array[ i + end ][ 'A' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'C' ) expr_array[ i + end ][ 'C' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'G' ) expr_array[ i + end ][ 'G' ] += read[ 'PPM' ];" << "\n";
        output << "                            if( read[ 'Tail' ].charAt(i) == 'T' ) expr_array[ i + end ][ 'T' ] += read[ 'PPM' ];" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( expr_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( expr_max < expr_array[ idx ][ 'PPM' ] ) expr_max = expr_array[ idx ][ 'PPM' ];" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( expr_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        expr_array[ idx ][ 'PPM' ] = expr_array[ idx ][ 'PPM' ] / expr_max * 100;" << "\n";
        output << "" << "\n";
        output << "                        var total = expr_array[ idx ][ 'P' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'A' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'C' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'G' ]" << "\n";
        output << "                                  + expr_array[ idx ][ 'T' ];" << "\n";
        output << "" << "\n";
        output << "                        if( total != 0 )" << "\n";
        output << "                        {" << "\n";
        output << "                            expr_array[ idx ][ 'P' ] = expr_array[ idx ][ 'P' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'A' ] = expr_array[ idx ][ 'A' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'C' ] = expr_array[ idx ][ 'C' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'G' ] = expr_array[ idx ][ 'G' ] / total;" << "\n";
        output << "                            expr_array[ idx ][ 'T' ] = expr_array[ idx ][ 'T' ] / total;" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( heat_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        if( heat_max < heat_array[ idx ] )" << "\n";
        output << "                            heat_max = heat_array[ idx ];" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $.each( heat_array, function( idx )" << "\n";
        output << "                    {" << "\n";
        output << "                        heat_array[ idx ] = heat_array[ idx ] / heat_max * 100;" << "\n";
        output << "                        if( heat_array[ idx ] == 0 ) heat_array[ idx ] = 0.01;" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $( 'body' ).append( \\\"<div id='chart'></div>\\\" );" << "\n";
        output << "                    $( '#chart' ).css({" << "\n";
        output << "                        'margin-top': '10px'," << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
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
        output << "                        });" << "\n";
        output << "" << "\n";
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
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"T'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"G'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"C'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"A'></div>\\\" );" << "\n";
        output << "                        $( '#expr' + idx ).append( \\\"<div id='expr\\\" + idx + \\\"P'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                        $( '#expr' + idx + 'T' ).css({ 'height': ( expr[ 'PPM' ] * expr[ 'T' ]) + 'px', 'background': colorT });" << "\n";
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
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabeltop'>-\\\" + expr_max.toFixed(2) + '</div>' );" << "\n";
        output << "                    $( '#exprlabeltop' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': 0 + shift_top + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabelmid'>-\\\" + ( expr_max /2 ).toFixed(2) + '</div>' );" << "\n";
        output << "                    $( '#exprlabelmid' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': 50 + shift_top + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#exprlabel' ).append( \\\"<div id='exprlabeldown'>-0.00</div>\\\" );" << "\n";
        output << "                    $( '#exprlabeldown' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': 100 + shift_top + 'px'," << "\n";
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
        output << "                        'pointer-events': 'none'," << "\n";
        output << "                        'height': '25px'," << "\n";
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
        output << "                            $( '#seqs' ).append( \\\"<a id='seq\\\" + i + \\\"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\\\" + seeds[i][ 'Name' ] + \\\"'></a>\\\" );" << "\n";
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
        output << "                    }" << "\n";
        output << "" << "\n";
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
        output << "                        'top': seq_height + exp_height + shift_top + 'px'," << "\n";
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
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.25 ) + seq_height + exp_height + shift_top + 'px'," << "\n";
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
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.5 ) + seq_height + exp_height + shift_top + 'px'," << "\n";
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
        output << "                            'top': ( datas.filter( Boolean ).length * read_height * 0.75 ) + seq_height + exp_height + shift_top + 'px'," << "\n";
        output << "                            'float': 'left'" << "\n";
        output << "                            });" << "\n";
        output << "                    }" << "\n";
        output << "" << "\n";
        output << "                    $( '#readslabel' ).append( \\\"<div id='readslabeldown'>-0.00</div>\\\" );" << "\n";
        output << "                    $( '#readslabeldown' ).css({" << "\n";
        output << "                        'font-family': 'Source Code Pro, monospace'," << "\n";
        output << "                        'font-size': '15px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': ( datas.filter( Boolean ).length * read_height ) + seq_height + exp_height + shift_top + 'px'," << "\n";
        output << "                        'float': 'left'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
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
        output << "                            if( i in MDarray ) MDseed += ( i -1 ) + MDarray[i];" << "\n";
        output << "" << "\n";
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
        output << "                            $( '#selectedseed' + last_idx ).append( \\\"<a id='read\\\" + last_idx + \\\"linked' target='_blank' href='../LenDist/index.php?TSV_File=$TSV_File_Name&annotation_select=\\\"" << "\n";
        output << "                                + seeds[ last_idx + 4 ][ 'Name' ]" << "\n";
        output << "                                + ( MDseed == '' ? '' : ( '|' + MDseed ))" << "\n";
        output << "                                + ( read[ 'isRMSK' ] == 'N' ? '' : '!' )" << "\n";
        output << "                                + \\\"'></a>\\\" );" << "\n";
        output << "" << "\n";
        output << "                            $( '#read' + last_idx + 'linked' ).css({ 'text-decoration': 'none' });" << "\n";
        output << "                        }" << "\n";
        output << "" << "\n";
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
        output << "                            if( read[ 'Tail' ].charAt(i) == 'T' ) $( '#read' + idx + 'pm' + i ).css({ 'color': colorT });" << "\n";
        output << "" << "\n";
        output << "                            if( i != read[ 'Tail' ].length -1 )" << "\n";
        output << "                            {" << "\n";
        output << "                                $( '#read' + idx + 'pm' + i ).css({" << "\n";
        output << "                                    'width': spc_num + 'px'," << "\n";
        output << "                                    'float': 'left'" << "\n";
        output << "                                    });" << "\n";
        output << "                            }" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "" << "\n";
        output << "                    $( '#chart' ).css({" << "\n";
        output << "                        'width': (( json[ 'Sequence' ].length + extra + lb_width )* spc_num ) + 'px'," << "\n";
        output << "                        'height': reads_count * read_height + seq_height + exp_height + 'px'" << "\n";
        output << "                        });" << "\n";
        output << "" << "\n";
        output << "                    $( '#chart' ).append( \\\"<div id='line'></div>\\\" );" << "\n";
        output << "" << "\n";
        output << "                    $( '#line' ).css({" << "\n";
        output << "                        'height': $( '#chart' ).height() + 2 + 'px'," << "\n";
        output << "                        'width': '1px'," << "\n";
        output << "                        'position': 'absolute'," << "\n";
        output << "                        'top': '40px'," << "\n";
        output << "                        'left': '-1px'," << "\n";
        output << "                        'background-color': 'indianred'," << "\n";
        output << "                        'pointer-events': 'none'" << "\n";
        output << "                        });" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
        output << "                var selectedread = $( '#chart' );" << "\n";
        output << "                var selectedseed = $( '#chart' );" << "\n";
        output << "" << "\n";
        output << "                $( document ).mousemove( function( event )" << "\n";
        output << "                {" << "\n";
        output << "                    $( '#line' ).css({" << "\n";
        output << "                        'left': ( event.pageX > $( '#chart' ).width() ? -1 : event.pageX ) + 'px'" << "\n";
        output << "                        });" << "\n";
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
        output << "                            if( selectedread.attr( 'title' ).substr( 0, 4 ) == 'expr' )" << "\n";
        output << "                            {" << "\n";
        output << "                                data = exprs[ selectedread.attr( 'title' ).substr( 4 ) ];" << "\n";
        output << "                                tip += 'PPM: ' + ( data[ 'PPM' ] * expr_max / 100 ).toFixed(2);" << "\n";
        output << "" << "\n";
        output << "                                if( data[ 'A' ] > 0.001 ) tip += '</br>Tail A: ' + ( data[ 'A' ] * 100 ).toFixed(2) + '%';" << "\n";
        output << "                                if( data[ 'C' ] > 0.001 ) tip += '</br>Tail C: ' + ( data[ 'C' ] * 100 ).toFixed(2) + '%';" << "\n";
        output << "                                if( data[ 'G' ] > 0.001 ) tip += '</br>Tail G: ' + ( data[ 'G' ] * 100 ).toFixed(2) + '%';" << "\n";
        output << "                                if( data[ 'T' ] > 0.001 ) tip += '</br>Tail T: ' + ( data[ 'T' ] * 100 ).toFixed(2) + '%';" << "\n";
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
        output << "                                tip += 'Annotation: ' + anno;" << "\n";
        output << "                                tip += \".( $Data_Array[0][5] != '.'? ( \"'-' + data[ 'Arm' ]\" ) : \"''\" ).\";" << "\n";
        output << "" << "\n";
        output << "                                tip += '_' + sequ.substr(( data[ 'Index' ] + 1 ), 7 );" << "\n";
        output << "                                tip += ( MDseed == '' ? '' : ( '|<b style=\\\"color:red;\\\">' + MDseed + '</b>' ));" << "\n";
        output << "" << "\n";
        output << "                                tip += ( data[ 'isRMSK' ] == 'N' ? '' : ' (RMSK)' ) + '</br>';" << "\n";
        output << "                                tip += 'Sequence: ';" << "\n";
        output << "" << "\n";
        output << "                                for( var i = 0; i < seq.length; ++i )" << "\n";
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
        output << "                                    selectedseed.css({" << "\n";
        output << "                                        'border': 'none'," << "\n";
        output << "                                        });" << "\n";
        output << "                                }" << "\n";
        output << "                            );" << "\n";
        output << "                            $( '#selectbox' ).css({ 'display': 'none' });" << "\n";
        output << "                            selectedread.css({ 'border': 'none' });" << "\n";
        output << "                            selectedseed.css({" << "\n";
        output << "                                'border': 'none'," << "\n";
        output << "                                });" << "\n";
        output << "                        }" << "\n";
        output << "                    });" << "\n";
        output << "                });" << "\n";
        output << "" << "\n";
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
