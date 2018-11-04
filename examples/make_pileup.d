// BioD depends on stream.d which is no longer included with phobos.
// To run this example from this directory,
// Clone the undead repository with
// git clone https://github.com:dlang/undeaD.git at the appropriate location and ensure
// it is available on your path
// Run example: rdmd -I.. -I../location_of_undead/src make_pileup.d

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.read : compareCoordinates;

import std.range;
import std.algorithm;
import std.datetime;
import std.stdio;
import std.array;

void main() {

  auto bam = new BamReader("../test/data/illu_20_chunk.bam");

  auto pileup = makePileup(bam.reads, true);

  // Reads in every pileup column are sorted by coordinate.
  // Therefore the following holds:
  assert(equal(
        pileup.map!(column => column.reads.equalRange(column.position))()
        .joiner(),
        bam.reads));

  // There is also easier and faster way to get reads starting at the column:
//  pileup = makePileup(bam.reads, true); // (initialize pileup engine again)
 // assert(equal(
 //       pileup.map!(column => column.reads_starting_here)()
  //      .joiner(),
  //      bam.reads));
}
