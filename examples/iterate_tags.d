// BioD depends on stream.d which is no longer included with phobos.
// To run this example from this directory,
// Clone the undead repository with
// git clone https://  /undead/undead.git at the appropriate location and ensure
// it is available on your path
// Run example: rdmd -I.. -I../location_of_undead/src iterate_tags.d

import bio.bam.reader, bio.bam.tagvalue;
import std.stdio, std.datetime;
import std.conv : to;

void main() {

  auto bam = new BamReader("../test/data/b7_295_chunk.bam");

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
  auto v3 = to!long(v); auto v4 = to!string(v); auto v5 = to!float(v);

  // With strings and arrays there is an unsafe but faster way...
  v = read["FZ"];
  StopWatch sw;

  // even with -O -release -inline this is slow
  sw.start; auto fz1 = to!(ushort[])(v); sw.stop();
  writeln("  safe conversion: ", sw.peek().usecs, "μs");
  sw.reset();

  // this works because v starts in memory with a union
  sw.start(); auto fz2 = *(cast(ushort[]*)(&v)); sw.stop();
  writeln("unsafe conversion: ", sw.peek().usecs, "μs");

  assert(fz1 == fz2);
}
