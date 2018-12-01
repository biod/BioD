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
module bio.std.file.fasta;

import std.file;
import std.exception;
import std.algorithm;
import std.string;

struct FastaRecord {
    string header;
    string sequence;
}

auto fastaRecords(string filename) {

    static auto toFastaRecord(S)(S str) {
        auto res = findSplit(str, "\n");
        auto header = res[0];
        auto seq = res[2];
        return FastaRecord(header, removechars(seq, "\n"));
    }

    string text = cast(string)std.file.read(filename);

    enforce(text.length > 0 && text[0] == '>');
    text = text[1 .. $];

    auto records = splitter(text, '>');
    return map!toFastaRecord(records);
}
