
import bio.std.hts.bam.reader;
import bio.std.hts.bam.pileup;
import bio.std.hts.bam.read : compareCoordinates;

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
 pileup = makePileup(bam.reads, true); // (initialize pileup engine again)
 assert(equal(
       pileup.map!(column => column.reads_starting_here)()
      .joiner(),
      bam.reads));
}
