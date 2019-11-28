
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
import core.memory;

void main() {
  BamReader bam;
  auto fn = buildPath(dirName(__FILE__), "data", "ex1_header.bam");
  foreach (i; 0 .. 60)
  {
    bam = new BamReader(fn);
    assert(bam.header.format_version == "1.3");
  }
  // Time to kick in GC
  // stderr.writeln("**** Calling GC final");
  // GC.collect();
  // stderr.writeln("**** Past calling GC final");
}
