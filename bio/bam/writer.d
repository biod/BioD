/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.bam.writer;

import bio.bam.referenceinfo;
import bio.sam.header;
import bio.bam.constants;
import bio.core.bgzf.outputstream;
import bio.core.utils.stream;
import bio.core.utils.switchendianness;

import std.parallelism;
import std.exception;
import std.stream;
import std.traits;
import std.system;

/** Class for outputting BAM.
    $(BR)
    Compresses BGZF blocks in parallel.
    Tries to write reads so that they don't cross BGZF block borders.
    $(BR)
    Usage is very simple, see example below.

    Example:
    --------------------------------------
    import bio.bam.writer, bio.bam.reader;
    ...
    auto src = new BamReader("in.bam");
    auto dst = new BamWriter("out.bam", 9); // maximal compression
    scope (exit) dst.finish();              // close the stream at exit
    dst.writeSamHeader(src.header);         // copy header and reference sequence info
    dst.writeReferenceSequenceInfo(src.reference_sequences);
    foreach (read; src.reads) {
        if (read.mapping_quality > 10)      // skip low-quality reads
            dst.writeRecord(read);
    }
    --------------------------------------
    */
final class BamWriter {

    /// Creates new BAM writer outputting to file or $(I stream).
    /// Automatically writes BAM magic number (4 bytes).
    ///
    /// Params:
    ///     compression_level  = compression level, must be in range -1..9
    ///     task_pool          = task pool to use for parallel compression
    ///     buffer_size        = size of BgzfOutputStream buffer
    this(std.stream.Stream stream, 
         int compression_level=-1,
         std.parallelism.TaskPool task_pool=std.parallelism.taskPool,
         size_t buffer_size=0) 
    {
        _stream = new BgzfOutputStream(stream, compression_level, 
                                       task_pool, buffer_size);
        writeString(BAM_MAGIC);
    }

    /// ditto
    this(string filename,
         int compression_level=-1,
         std.parallelism.TaskPool task_pool=std.parallelism.taskPool)
    {
        auto filestream = new bio.core.utils.stream.File(filename, "wb+");
        this(filestream, compression_level, task_pool);
    }

    package void writeByteArray(const(ubyte[]) array) {
        _stream.writeExact(array.ptr, array.length);
    }

    package void writeString(string str) {
        writeByteArray(cast(ubyte[])str);
    }

    package void writeInteger(T)(T integer) if (isIntegral!T)
    {
        T num = integer;
        static if (T.sizeof != 1) {
            if (std.system.endian != Endian.littleEndian) {
                switchEndianness(&num, T.sizeof);
            }
        }
        _stream.writeExact(&num, T.sizeof);
    }

    /// Writes SAM header. Should be called after construction.
    void writeSamHeader(bio.sam.header.SamHeader header) {
        auto text = header.text;
        writeInteger(cast(int)text.length);
        writeString(text);
    }

    /// Writes reference sequence information. Should be called after
    /// dumping SAM header. Writer will store this array to use later for
    /// resolving read reference IDs to names.
    ///
    /// Flushes current BGZF block.
    void writeReferenceSequenceInfo(const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences)
    {
        _reference_sequences = reference_sequences;

        writeInteger(cast(int)reference_sequences.length);
        foreach (sequence; reference_sequences) {
            writeInteger(cast(int)(sequence.name.length + 1));
            writeString(sequence.name);
            writeInteger(cast(ubyte)'\0');
            writeInteger(cast(int)sequence.length);
        }

        _stream.flushCurrentBlock();
    }

    /// Writes BAM read. Throws exception if read reference ID is out of range.
    void writeRecord(R)(R read) {
        enforce(read.ref_id == -1 || read.ref_id < _reference_sequences.length,
                "Read reference ID is out of range");

        auto read_size = read.size_in_bytes;
        if (read_size + _current_size > BGZF_BLOCK_SIZE) {
            _stream.flushCurrentBlock();
            read.write(this);
            _current_size = read_size;
        } else {
            read.write(this);
            _current_size += read_size;
        }
    }

    /// Flushes current BGZF block.
    void flush() {
        _stream.flush();
    }

    /// Flushes buffer and closes output stream. Adds BAM EOF block automatically.
    void finish() {
        _stream.close();
    }

    private {
        BgzfOutputStream _stream;
        const(ReferenceSequenceInfo)[] _reference_sequences;
        size_t _current_size; // number of bytes written to the current BGZF block
    }
}
