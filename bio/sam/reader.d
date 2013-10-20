/*
    This file is part of BioD.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
module bio.sam.reader;

import bio.bam.abstractreader;
import bio.sam.header;
import bio.bam.read;
import bio.bam.referenceinfo;
import bio.core.utils.outbuffer;
import bio.core.utils.range;

import bio.core.utils.bylinefast;
alias ByLineFast _LineRange;

version(DigitalMars) {
    import bio.sam.utils.recordparser;
} else {
    import bio.sam.utils.fastrecordparser;
}

import std.stdio;
import std.array;
import std.string;
import std.range;
import std.algorithm;
import std.typecons;
import std.parallelism;
import std.c.string;

BamRead _parseSamRecord(Tuple!(char[], SamReader, OutBuffer) t) {
    auto r = parseAlignmentLine(cast(string)t[0], t[1]._header, t[2]);
    auto storage = uninitializedArray!(ubyte[])(r.raw_data.length);
    storage[] = r.raw_data[];
    BamRead result;
    result.raw_data = storage;
    result.associateWithReader(t[1]);
    return result;
}

BamRead[] _parseSamRecords(Tuple!(string[], SamReader) t) {
    auto lines = t[0];
    auto reader = t[1];

    BamRead[] reads;
    reads.reserve(lines.length);

    enum chunk_size = 65536;
    auto b = new OutBuffer(chunk_size);
    ubyte[] buffer = uninitializedArray!(ubyte[])(chunk_size);
    size_t used;
    
    foreach (line; lines) {
        auto r = parseAlignmentLine(line, reader.header, b);
        auto len = r.raw_data.length;
        if (len + used > buffer.length) {
            used = 0;
            buffer = uninitializedArray!(ubyte[])(max(len, chunk_size));
        }
        buffer[used .. used + len] = r.raw_data[];
        BamRead read;
        ubyte[] raw = buffer[used .. used + len];
        read.raw_data = raw;
        used += len;
        read.associateWithReader(reader);
        reads ~= read;
    }
    return reads;
}

private {
    extern(C) size_t lseek(int, size_t, int);
    bool isSeekable(ref File file) {
        return lseek(file.fileno(), 0, 0) != ~0;
    }

    struct SamChunks {
        private {
            _LineRange _lines;
            size_t _max_sz;
            ubyte[] _data;
            size_t _used;
            string[] _front;
            bool _empty;

            void getNextChunk() {
                if (_lines.empty) {
                    _empty = true;
                    return;
                }

                assert(_used == 0);
                auto result = Appender!(string[])();
                result.reserve(_max_sz / 128);

                while (!_lines.empty) {
                    auto _line = _lines.front;
                    auto len = _line.length;
                    if (_line.length + _used > _data.length) {
                        _data = uninitializedArray!(ubyte[])(max(_line.length, _max_sz));
                        _used = 0;
                        break; // this line goes to the next chunk
                    }

                    std.c.string.memcpy(_data.ptr + _used, _line.ptr, _line.length);
                    result.put(cast(string)(_data[_used .. _used + len]));
                    _used += len;
                               
                    _lines.popFront();
                }

                _front = result.data;
            }
        }

        bool empty() @property { return _empty; }
        string[] front() @property { return _front; }
        void popFront() { getNextChunk(); }

        this(_LineRange lines, size_t max_size) {
            assert(max_size > 2^^16);
            _lines = lines;
            _max_sz = max_size;
            if (_lines.empty) {
                _empty = true;
            } else {
                _data = uninitializedArray!(ubyte[])(max(_lines.front.length, _max_sz));
                getNextChunk();
            }
        }
    }
}

///
class SamReader : IBamSamReader {

    ///
    this(string filename) {
        _file = File(filename);
        _filename = filename;
        _seekable = _file.isSeekable();
        _initializeStream();
    }

    ///
    bio.sam.header.SamHeader header() @property {
        return _header;
    }

    ///
    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const {
        return _reference_sequences;
    }

    /// Reads in SAM file.
    auto reads() @property {
        
        _LineRange lines = _lines;
        if (_seekable) {
            if (_filename !is null) {
                File file = File(_filename);
                lines = ByLineFast(file);
            } else {
                _file.seek(0);
                lines = ByLineFast(_file);
            }
            auto dummy = lines.front;
            for (int i = 0; i < _lines_to_skip; i++)
                lines.popFront();
        }

        auto b = new OutBuffer(262144);
        return lines.zip(repeat(this), repeat(b)).map!_parseSamRecord();
        //   return SamChunks(lines, 1_024_576).zip(repeat(this))
        //    .map!_parseSamRecords().joiner();
    }

    ///
    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads);
    }

    /// Filename
    string filename() @property const {
        return _filename;
    }
private:

    File _file;
    bool _seekable;
    string _filename;
    _LineRange _lines;
    ulong _lines_to_skip;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream() {
        auto header = Appender!(char[])(); 

        _lines = ByLineFast(_file);

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
                header.put('\n');
                _lines_to_skip += 1;
                _lines.popFront();
            } else {
                break;
            }
        }

        _header = new SamHeader(cast(string)(header.data));

        _reference_sequences = new ReferenceSequenceInfo[_header.sequences.length];
        foreach (sq; _header.sequences) {
            auto seq = ReferenceSequenceInfo(sq.name, sq.length);
            _reference_sequences[_header.getSequenceIndex(seq.name)] = seq;
        }
    }
}
