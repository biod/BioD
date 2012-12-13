// run example: rdmd -I.. example5.d

import bio.bam.reader;
import bio.bam.read : compareCoordinates;
import bio.bam.pileup;

import std.algorithm, std.conv, std.stdio;

void main() {
    auto bam1 = new BamReader("../test/data/illu_20_chunk.bam");
    auto bam2 = new BamReader("../test/data/ion_20_chunk.bam");

    // this is how to iterate over several datasets simultaneously
    // (internally, nWayUnion maintains a binary heap)
    auto reads = nWayUnion!compareCoordinates([bam1.reads, bam2.reads]);

    auto pileup = makePileup(reads, // ANY range of reads is acceptable
                             true,  // use MD tags
                             32_000_083,
                             32_000_089);

    foreach (column; pileup) {
        writeln("Column position: ", column.position);
        writeln("    Ref.base: ", column.reference_base); // extracted from MD tags
        writeln("    Coverage: ", column.coverage);

        writeln("    ", column.reads // bases from Illumina dataset
                              .filter!(read => to!string(read["RG"]).startsWith("ERR"))()
                              .map!(read => read.current_base)(),

                "    ", column.reads // bases from IonTorrent dataset
                              .filter!(read => to!string(read["RG"]).startsWith("66A0Q"))()
                              .map!(read => read.current_base)());
    }
}
