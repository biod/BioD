/* Run this example from this directory by typing this from the command line:

// run with rdmd
rdmd -I.. iterate_tags.d

// A debug compile with dmd 
dmd -i -I.. iterate_tags.d

// optimised compile 
dmd -i -O -release -inline -boundscheck=off -I.. iterate_tags.d

*/


import bio.std.hts.bam.reader;
import bio.std.hts.bam.tagvalue;
import std.stdio;
import std.datetime.stopwatch : benchmark,StopWatch;
import std.conv : to;

void main() {

    // add a file handling exception in the example 

    auto bam = new BamReader("../test/data/b7_295_chunk.bam");
    //  auto bam = new BamReader("/Users/george/HH_bam_files/H_3801_01_03.final.bam");

    auto read = bam.reads.front; // take first read

    // iterating all tags
    foreach (tag, value; read)
        writeln(tag, ": ", value);

    // taking value of tag
    Value v = read["XS"];

    // Usually, it will be converted to some type right away.
    auto v2 = to!int(v);

    // It is not necessary to know exact value type as in BAM.
    // If it can be converted to specified type, that's fine.
    // Otherwise, an exception will be thrown.
    auto v3 = to!long(v); 
    auto v4 = to!string(v); 
    auto v5 = to!float(v);

    // With strings and arrays there is an unsafe but faster way...
    v = read["FZ"];

    StopWatch sw;

    // even with -O -release -inline this is slow
    sw.start; 
    auto fz1 = to!(ushort[])(v);
    sw.stop();

    writeln("  safe conversion: ", sw.peek.total!"usecs", "μs");
    sw.reset();

    // this works because v starts in memory with a union
    sw.start();
    auto fz2 = *(cast(ushort[]*)(&v)); 
    sw.stop();

    writeln("unsafe conversion: ", sw.peek.total!"usecs", "μs");

    assert(fz1 == fz2);
}
