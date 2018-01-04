// BioD depends on stream.d which is no longer included with phobos.
// To run this example from this directory,
// Clone the undead repository with
// git clone https://  /undead/undead.git at the appropriate location and ensure
// it is available on your path
// Run example: rdmd -I.. -I../location_of_undead/src transverse_multiple_bam_files.d

import bio.bam.multireader;
import bio.bam.read : compareCoordinates;
import bio.bam.pileup;

import std.algorithm, std.conv, std.stdio;

void main() {
  // multiple BAM files can be traversed simultaneously (if they can be merged)
  auto bam = new MultiBamReader(["../test/data/illu_20_chunk.bam",
      "../test/data/ion_20_chunk.bam"]);

  auto pileup = makePileup(bam.reads, // ANY range of reads is acceptable
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
