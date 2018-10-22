// BioD depends on stream.d which is no longer included with phobos.
// To run this example from this directory,
// Clone the undead repository with
// git clone https://github.com:dlang/undeaD.git at the appropriate location and ensure
// it is available on your path
// Run example: rdmd -I.. -I../location_of_undead/src create_bam_from_scratch.d

// this example shows how to create BAM files from scratch
import bio.bam.read;
import bio.bam.referenceinfo;
import bio.sam.header;
import bio.bam.reader;
import bio.bam.writer;
import undead.stream;
import std.stdio;
import bio.bam.cigar;

void main() {
    auto header = new SamHeader(); // First, create SAM header 
    RgLine rg;                     // and fill it with information.
    rg.identifier = "RG007";       
    rg.platform = "ILLUMINA";      // Of course, you can modify header
    header.read_groups.add(rg);    // provided by BamReader object.

    auto reference = ReferenceSequenceInfo("contig123", 4321);
    auto stream = new MemoryStream();
    auto writer = new BamWriter(stream);
    writer.writeSamHeader(header);                  // start writing BAM from header
    writer.writeReferenceSequenceInfo([reference]); // and reference sequence info

    auto read = BamRead("readName001", "ACTGATGAAC",
                        [CigarOperation(7, 'M'), CigarOperation(1, 'I'), CigarOperation(2, 'S')]);
    read.base_qualities = [38, 34, 33, 35, 28, 39, 25, 19, 21, 17];
    read.ref_id = 0;
    read.mapping_quality = 46;
    read.is_unmapped = false;
    read["RG"] = "RG007";
    // ... set many more fields, flags, and tags
    read.position = 2345; // 0-based, in SAM output you will see 2346
    read.strand = '-'; // same as read.is_reverse_strand = true
    writer.writeRecord(read); // BGZF blocks are seamlessly compressed in parallel
    writer.flush(); // in practice, one would call finish() method
    stream.seekSet(0); // but here we will read from the stream

    auto reader = new BamReader(stream);
    write(reader.header.text); // serialized header already contains newline
    writefln("%j", reader.reads.front); // prints record in JSON format
}
