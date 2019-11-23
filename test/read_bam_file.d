import bio.std.hts.bam.reader;
import bio.std.hts.bam.pileup;
import std.stdio;

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

    auto fn = buildPath(dirName(__FILE__), "data", "tags.bam");
    auto bf = new BamReader(fn);
    foreach (read; bf.reads) {
        auto line = read.to!string();
        auto read2 = parseAlignmentLine(line, bf.header);
        if (read != read2 && isValid(read)) {
            stderr.writeln(read.name);
        }
        assert(read == read2 || !isValid(read));
    }

    fn = buildPath(dirName(__FILE__), "data", "ex1_header.bam");
    bf = new BamReader(fn);
    auto readx = bf.reads;
    readx.front;
    readx.popFront();
    readx.popFront();

        {
      fn = buildPath(dirName(__FILE__), "data", "b.sam");
      auto sam = new SamReader(fn);
      auto writer = new BamWriter("/dev/null", 0);
      writer.writeSamHeader(sam.header);
      writer.writeReferenceSequenceInfo(sam.reference_sequences);
      foreach (r; sam.reads)
        writer.writeRecord(r);
      writer.finish();
    }
    {
      fn = buildPath(dirName(__FILE__), "data", "bins.bam");
      bf = new BamReader(fn);

      void compareWithNaiveApproach(int beg, int end) {

        auto refseq = array(bf["large"][beg .. end]);

        auto naive = array(filter!((BamRead a) {
                         return a.ref_id != -1 &&
                                bf.reference(a.ref_id).name == "large" &&
                                a.position < end &&
                                a.position + a.basesCovered() > beg; })
                            (bf.reads!withoutOffsets));
        if (!equal(naive, refseq)) {
            stderr.writeln(beg);
            stderr.writeln(end);
            stderr.writeln(array(map!"a.name"(refseq)));
            stderr.writeln(array(map!"a.name"(naive)));
        }
        assert(equal(refseq, naive));
      }


      // Time to kick in GC
      import core.memory;
      stderr.writeln("**** Calling GC");
      GC.collect();
      stderr.writeln("**** Past calling GC");
    }
CigarOperation[] cigarFromString(string cigar) {
    return match(cigar, regex(`(\d+)([A-Z=])`, "g")).map!(m => CigarOperation(m[1].to!uint, m[2].to!char)).array;
}

{
    stderr.writeln("Running unittests...");
    // stderr.writeln("Testing extracting SAM header...");

    fn = buildPath(dirName(__FILE__), "data", "ex1_header.bam");
    bf = new BamReader(fn);
    assert(bf.header.format_version == "1.3");
    assert(bf.header.sorting_order == SortingOrder.coordinate);
    assert(bf.header.sequences.length == 2);
    assert(bf.header.getSequenceIndex("chr1") == 0);
    assert(bf.header.sequences["chr2"].length == 1584);

    fn = buildPath(dirName(__FILE__), "data", "bins.bam");
    bf = new BamReader(fn);
    assert(bf.header.sorting_order == SortingOrder.unknown);
    assert(bf.header.sequences.length == 3);
    assert(bf.header.read_groups.length == 0);
    assert(bf.header.getSequenceIndex("large") == 2);
    assert(bf.header.sequences["small"].length == 65536);

    // Time to kick in GC
    import core.memory;
    stderr.writeln("**** Calling GC");
    GC.collect();
    stderr.writeln("**** Past calling GC");

    {
    // stderr.writeln("Testing alignment parsing...");
    fn = buildPath(dirName(__FILE__), "data", "ex1_header.bam");
    bf = new BamReader(fn);
    auto reads2 = bf.reads;
    auto read = reads2.front;
    assert(equal(read.sequence, "CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA"));
    assert(equal(map!"cast(char)(a + 33)"(read.base_qualities),
                "<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;"));
    assert(bf.reference(read.ref_id).name == "chr1");
    assert(read.name == "EAS56_57:6:190:289:82");
    assert(read.flag == 69);
    assert(read.position == 99);
    assert(read.mapping_quality == 0);
    reads2.popFront();
    reads2.popFront();
    assert(reads2.front.cigarString() == "35M");
    assert(reads2.front.to!string() == "EAS51_64:3:190:727:308	99	chr1	103	99	35M	=	263	195	GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG	<<<<<<<<<<<<<<<<<<<<<<<<<<<::<<<844	MF:i:18	Aq:i:73	NM:i:0	UQ:i:0	H0:i:1	H1:i:0");
    assert(bf.header.getSequenceIndex("chr1") == read.ref_id);
    }

    assert(bf.reads.front.name == "EAS56_57:6:190:289:82");
    stderr.writeln("**** Calling GCx");
    GC.collect();
    stderr.writeln("**** Past calling GCx");

    // stderr.writeln("Testing tag parsing...");
    fn = buildPath(dirName(__FILE__), "data", "tags.bam");
    bf = new BamReader(fn);
    foreach (alignment; bf.reads) {
        auto name = alignment.name;
        assert(name[0..4] == "tag_");
        char[] tag;
        name = name[4..$];
        while (name[0] != ':') {
            tag ~= name[0];
            name = name[1..$];
        }
        name = name[1..$];
        auto value = alignment[tag.idup].toSam();
        if (name != value) {
            stderr.writeln("tag: ", tag, "\tname: ", name, "\tvalue: ", value);
            stderr.writeln("value bam_typeid: ", alignment[tag.idup].bam_typeid);
        }

        assert(name == value);
    }

    // stderr.writeln("Testing exception handling...");
    fn = buildPath(dirName(__FILE__), "data", "duplicated_block_size.bam");
    assertThrown!BgzfException(new BamReader(fn));
    fn = buildPath(dirName(__FILE__), "data", "no_block_size.bam");
    assertThrown!BgzfException(new BamReader(fn));
    fn = buildPath(dirName(__FILE__), "data", "wrong_extra_gzip_length.bam");
    assertThrown!BgzfException(new BamReader(fn));
    fn = buildPath(dirName(__FILE__), "data", "wrong_bc_subfield_length.bam");
    assertThrown!BgzfException(reduce!"a+b.sequence_length"(0, (new BamReader(fn)).reads!withoutOffsets));
    fn = buildPath(dirName(__FILE__), "data", "corrupted_zlib_archive.bam");
    import bio.core.utils.zlib;
    assertThrown!ZlibException(walkLength((new BamReader(fn)).reads));
    // stderr.writeln("Testing random access...");


  }
}
