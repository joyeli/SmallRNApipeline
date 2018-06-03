#include <iostream>
#include <AGO/format/md_sam.hpp>


int main( int argc, char** argv )
{
    std::string samstr;
    samstr = "test\t0\tchr14\t43718720\t251\t20M4S\t*\t0\t0\tAGCGGAGTGTATAATGGCTTAAAA\tIIIIIIIIIIIII*III%IIII%1\tNH:i:1\tMD:Z:7A9T2\tTL:Z:AAAA";
    // samstr = "test\t0\tchr14\t43718720\t251\t20M4S\t*\t0\t0\tAGCGGAGTGTATAATGGCTTAAAA\tIIIIIIIIIIIII*III%IIII%1\tNH:i:1\tMD:Z:7A9T2";
    // samstr = "test\t0\tchr14\t43718720\t251\t20M4S\t*\t0\t0\tAGCGGAGTGTATAATGGCTTAAAA\tIIIIIIIIIIIII*III%IIII%1\tNH:i:1\tTL:Z:AAAA";
    // samstr = "test\t0\tchr14\t43718720\t251\t20M4S\t*\t0\t0\tAGCGGAGTGTATAATGGCTTAAAA\tIIIIIIIIIIIII*III%IIII%1\tNH:i:1";

    ago::format::MDSam<> md_sam( samstr );

    return 0;
}
