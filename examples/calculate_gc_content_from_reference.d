/*
   To run this example from this directory:
   rdmd -I.. calculate_gc_content_from_reference.d

// compile this example with

// debug version
dmd -i -I.. calculate_gc_content_from_reference.d

// optimised version
dmd -i -O -release -inline -boundscheck=off -I.. calculate_gc_content_from_reference.d

*/
import bio.std.hts.bam.reader;
import bio.std.hts.bam.md.reconstruct : dna;

import std.stdio;
import std.datetime.stopwatch: benchmark, StopWatch;
import std.range;
import std.array;

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

    StopWatch sw; // for a range of reads, minimum number of MD tags is parsed

    sw.start(); 
    walkLength(dna(reads));
    sw.stop();

    writeln("extracting reference from range of reads: ", sw.peek.total!"usecs", "μs");

    sw.reset(); 
    sw.start(); 
    foreach (read; reads) walkLength(dna(read)); 
    sw.stop();

    writeln("extracting reference from each read:     ", sw.peek.total!"usecs", "μs");
}
