/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    BioD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    BioD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

/// Class for outputting BAM.
/// Tries to write reads so that they don't cross BGZF block borders.
final class BamWriter {

    /// Create new BAM writer outputting to $(D stream).
    /// Automatically writes BAM magic number (4 bytes).
    ///
    /// Params:
    ///     compression_level  = compression level, must be in range -1..9
    ///     task_pool          = task pool to use for parallel compression
    this(Stream stream, 
        int compression_level=-1,
        TaskPool task_pool=taskPool) 
    {
        _stream = new BgzfOutputStream(stream, compression_level, task_pool);
        writeString(BAM_MAGIC);
    }

    /// Create new BAM writer outputting to specified file.
    this(string filename,
         int compression_level=-1,
         TaskPool task_pool=taskPool)
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

    /// Write SAM header. Should be called after construction.
    void writeSamHeader(SamHeader header) {
        auto text = toSam(header);
        writeInteger(cast(int)text.length);
        writeString(text);
    }

    /// Write reference sequence information. Should be called after
    /// dumping SAM header. Writer will store this array to use later for
    /// resolving read reference IDs to names.
    ///
    /// Flushes current BGZF block.
    void writeReferenceSequenceInfo(ReferenceSequenceInfo[] reference_sequences)
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

    /// Write BAM read. Throws exception if read reference ID is out of range.
    void writeRecord(R)(R read) {
        enforce(read.ref_id < _reference_sequences.length,
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

    /// Flush data into stream.
    void flush() {
        _stream.flush();
    }

    /// Flush buffer and close output stream. Adds BAM EOF block automatically.
    void finish() {
        _stream.close();
    }

    private {
        BgzfOutputStream _stream;
        ReferenceSequenceInfo[] _reference_sequences;
        size_t _current_size; // number of bytes written to the current BGZF block
    }
}
