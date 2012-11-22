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
module bio.core.fasta;

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
