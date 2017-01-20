#include <iostream>
#include <gtest/gtest.h>
#include <pokemon/format/annotation_raw_bed.hpp>
#include <pokemon/annotator/annotation.hpp>
#include <pokemon/annotator/annotation_set.hpp>

using BedFileReaderImpl =  FileReader_impl<
      Bed
    , std::tuple< std::string, uint32_t, uint32_t, char, std::string, std::string >
    , SOURCE_TYPE::IFSTREAM_TYPE
>;

using AnnotationTrait = Annotation<
      BedFileReaderImpl
    , AnnoIgnoreStrand::NO_IGNORE
    , AnnoType::INTERSET
>;

using Annotations = AnnotationSet <
      AnnotationRawBed<>
    , AnnotationTrait
>;

TEST( annotator, error_test )
{
    std::vector< std::string > annotation_files{ "/home/joyel/AGO/bin/20170118/AnnoDB.bed" };
    Annotations annotator( annotation_files );
    AnnotationRawBed<> anno_rawbed;

    anno_rawbed.strand_ = '+';
    anno_rawbed.end_ = 59358338;
    anno_rawbed.chromosome_ = "chrY";
    anno_rawbed.annotation_info_ = {};
    anno_rawbed.is_filtered_ = 0;

    annotator.AnnotateAll( anno_rawbed );
}
