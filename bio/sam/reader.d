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

alias File.ByLine!(char, char) _LineRange;
BamRead _parseSamRecord(Tuple!(char[], SamReader, OutBuffer) t) {
    auto r = parseAlignmentLine(cast(string)t[0], t[1]._header, t[2]);
    r.associateWithReader(t[1]);
    return r;
}

private {
    extern(C) size_t lseek(int, size_t, int);
    bool isSeekable(ref File file) {
        return lseek(file.fileno(), 0, 0) != ~0;
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

    private alias File.ByLine!(char, char) LineRange;

    /// Reads in SAM file.
    auto reads() @property {
        
        LineRange lines = _lines;
        if (_seekable) {
            if (_filename !is null) {
                File file = File(_filename);
                lines = file.byLine();
            } else {
                _file.seek(0);
                lines = _file.byLine();
            }
            auto dummy = lines.front;
            for (int i = 0; i < _lines_to_skip; i++)
                lines.popFront();
        }

        auto buffer = new OutBuffer(8192);
        return lines.zip(repeat(this), repeat(buffer))
                    .map!_parseSamRecord().cached();
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
    LineRange _lines;
    ulong _lines_to_skip;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream() {
        auto header = Appender!(char[])(); 

        _lines = _file.byLine();

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
