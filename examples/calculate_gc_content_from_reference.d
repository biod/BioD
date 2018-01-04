// BioD depends on stream.d which is no longer included with phobos.
// To run this example from this directory,
// Clone the undead repository with
// git clone https://  /undead/undead.git at the appropriate location and ensure
// it is available on your path
// Run example: rdmd -I.. -I../location_of_undead/src calculate_gc_content_from_reference.d

import bio.bam.reader;
import bio.bam.md.reconstruct : dna;

import std.stdio, std.datetime, std.range, std.array;

void main() {
  auto bam = new BamReader("../test/data/b7_295_chunk.bam");

  // the sequence starts at first mapped base of first read
  auto reference = dna(bam.reads);

  int n_bases = 0, gc = 0;

  foreach (base; reference) {
    if (base == 'N') continue; // happens if coverage is zero
    if (base == 'G' || base == 'C') gc += 1;
    n_bases += 1;
  }

  writeln("total bases: ", n_bases);
  writeln("        GC%: ", cast(float)gc / n_bases);
  writeln("   sequence: ", reference);
  writeln("     #reads: ", walkLength(bam.reads));

  auto reads = array(bam.reads); // no I/O during measurements

  StopWatch sw; // for a range of reads, minumum number of MD tags is parsed
  sw.start(); walkLength(dna(reads)); sw.stop();
  writeln("extracting reference from range of reads: ", sw.peek().usecs, "μs");
  sw.reset(); 
  sw.start(); foreach (read; reads) walkLength(dna(read)); sw.stop();
  writeln("extracting reference from each read:     ", sw.peek().usecs, "μs");
}
