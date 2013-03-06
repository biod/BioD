/*
    This file is part of BioD.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>

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
module bio.maf.reader;

import bio.maf.block;
import bio.maf.parser;

import std.array;
import std.string;
import std.stdio;
import std.algorithm;

///
struct MafBlockRange {
    private {

        alias File.ByLine!(char, char) LineRange;
        File _f;
        LineRange _lines;

        bool _empty;
        MafBlock _front;

        void skipHeader() {
            if (!_lines.empty && _lines.front.startsWith("##maf"))
                _lines.popFront();
        }
    }

    this(string fn) {
        _f = File(fn);
        _lines = _f.byLine(KeepTerminator.yes);
        skipHeader();
        popFront();
    }

    ///
    bool empty() @property const {
        return _empty;
    }

    ///
    MafBlock front() @property {
        return _front;
    }

    ///
    void popFront() {
        auto block_data = Appender!(char[])();
        while (!_lines.empty && !_lines.front.chomp().empty) {
            block_data.put(_lines.front.dup);
            _lines.popFront();
        }
        if (block_data.data.empty) {
            _empty = true;
        } else {
            _front = parseMafBlock(cast(string)(block_data.data));
            if (!_lines.empty)
                _lines.popFront();
        }
    }
}


///
class MafReader {
 
    private string _fn;

    ///
    this(string filename) {
        _fn = filename;
    }

    ///
    string filename() @property const {
        return _fn;
    }

    ///
    MafBlockRange blocks() @property {
        return MafBlockRange(_fn);
    }
}
