
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.reader;
import bio.std.hts.bam.writer;
import bio.std.hts.sam.reader;
import bio.std.hts.sam.header;
import bio.std.hts.bam.md.core;
import bio.std.hts.bam.md.reconstruct;
import bio.std.hts.bam.pileup;
import bio.std.hts.bam.baseinfo;
import bio.std.hts.bam.validation.samheader;
import bio.std.hts.bam.validation.alignment;
import bio.std.hts.utils.samheadermerger;
import bio.std.hts.sam.utils.recordparser;
import bio.core.bgzf.block;
import bio.core.bgzf.inputstream;
import bio.core.bgzf.outputstream;
import bio.core.utils.roundbuf;
import bio.core.utils.range;
import bio.core.utils.tmpfile;
import bio.core.utils.stream;
import bio.core.sequence;
import bio.core.base;
import bio.core.tinymap;
import bio.core.utils.roundbuf;

import std.path;
import std.range;
import std.stdio;
// import undead.stream;
import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.math;
import std.typetuple;
import std.regex;

void main() {

    auto bam = new BamReader("test/data/ex1_header.bam");
    auto reads = bam["chr2"][150 .. 160]; // region chr2:149-158
    auto pileup = makePileup(reads,
            false,     // reads don't contain MD tags
            155, 158); // specify [start, end) interval

    foreach (column; pileup) {
        writeln("Reference position: ", column.position);
        writeln("    Coverage: ", column.coverage);
        writeln("    Reads:");

        foreach (read; column.reads) {
            writefln("%30s\t%s\t%.2d\t%s\t%2s/%2s\t%2s/%2s\t%10s\t%s %s",
                    read.name,
                    read.current_base,
                    read.current_base_quality,
                    read.cigar_operation,
                    read.cigar_operation_offset + 1, read.cigar_operation.length,
                    read.query_offset + 1, read.sequence.length,
                    read.cigarString(),
                    read.cigar_before, read.cigar_after);
        }
    }

 {
   auto fn = buildPath(dirName(__FILE__), "data", "b.sam");
   auto sam = new SamReader(fn);
   auto writer = new BamWriter("/dev/null", 0);
   writer.writeSamHeader(sam.header);
   writer.writeReferenceSequenceInfo(sam.reference_sequences);
   foreach (r; sam.reads)
     writer.writeRecord(r);
   writer.finish();
 }
 {
   auto fn = buildPath(dirName(__FILE__), "data", "bins.bam");
   auto bf = new BamReader(fn);
 }
   import core.memory;
   stderr.writeln("**** Calling GC");
   GC.collect();
   stderr.writeln("**** Past calling GC");
 {

   auto fn = buildPath(dirName(__FILE__), "data", "ex1_header.bam");
   auto bf = new BamReader(fn);
   assert(bf.header.format_version == "1.3");
   assert(bf.header.sorting_order == SortingOrder.coordinate);
   assert(bf.header.sequences.length == 2);
   assert(bf.header.getSequenceIndex("chr1") == 0);
   assert(bf.header.sequences["chr2"].length == 1584);

   auto fn2 = buildPath(dirName(__FILE__), "data", "bins.bam");
   auto bf2 = new BamReader(fn2);
   assert(bf2.header.sorting_order == SortingOrder.unknown);
   assert(bf2.header.sequences.length == 3);
   assert(bf2.header.read_groups.length == 0);
   assert(bf2.header.getSequenceIndex("large") == 2);
   assert(bf2.header.sequences["small"].length == 65536);
 }

 // Time to kick in GC
 import core.memory;
 stderr.writeln("**** Calling GC");
 GC.collect();
 stderr.writeln("**** Past calling GC");
}
