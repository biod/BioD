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
module bio.sam.reader;

import bio.bam.abstractreader;
import bio.sam.header;
import bio.bam.read;
import bio.bam.referenceinfo;

version(DigitalMars) {
    import bio.sam.utils.recordparser;
} else {
    import bio.sam.utils.fastrecordparser;
}

import std.stdio;
import std.array;
import std.string;

private {
    extern(C) size_t lseek(int, size_t, int);
}

///
class SamReader : IBamSamReader {

    ///
    this(string filename) {
        _file = File(filename);
        _filename = filename;
        _seekable = lseek(_file.fileno(), 0, 0) != ~0;
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

    static struct SamRecordRange {
        this(LineRange lines, ref SamHeader header, SamReader reader) {
            _header = header;
            _reader = reader;
            _line_range = lines;

            _build_storage = new AlignmentBuildStorage();
            _parseNextLine();
        }
        
        bool empty() @property {
            return _empty;
        }
        
        void popFront() @property {
            _line_range.popFront();
            _parseNextLine();
        }

        BamRead front() @property {
            return _current_alignment;
        }

        private {
            void _parseNextLine() {
                if (_line_range.empty) {
                    _empty = true;
                } else {
                    _current_alignment = parseAlignmentLine(cast(string)_line_range.front.dup,
                                                            _header,
                                                            _build_storage);
                    _current_alignment.associateWithReader(cast(IBamSamReader)_reader);
                }
            }

            LineRange _line_range;
            BamRead _current_alignment;
            bool _empty;
            SamHeader _header;
            SamReader _reader;
            AlignmentBuildStorage _build_storage;
        }
    }

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

        return SamRecordRange(lines, _header, this);
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
        auto header = appender!(char[])(); 

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
