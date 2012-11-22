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
module bio.bam.splitter;

import bio.bam.read;

import std.array;
import std.functional;
import std.range;
import std.traits;

/// Constructs range of chunks where sum { fn(read) | read in chunk }
/// does not exceed given number.
struct ReadRangeSplitter(R, alias fn) 
{
    this(R range, size_t threshold, bool split_by_ref) {
        _range = range;
        _size = threshold;
        _split_by_ref = split_by_ref;
        _appender = appender!(ElementType!R[])();
        getNextChunk();
    }

    private {
        R _range;
        bool _empty;
        bool _split_by_ref;
        size_t _size;
        Appender!(ElementType!R[]) _appender;
    }

    bool empty() @property {
        return _empty;
    }

    ElementType!R[] front() @property {
        return _appender.data.dup;
    }

    void popFront() {
        _appender.clear();
        getNextChunk(); 
    }

    private void getNextChunk() {
        if (_range.empty) {
            _empty = true;
            return;
        } 

        auto first_read = _range.front;
        _range.popFront();

        size_t total_size = first_read.size_in_bytes;
        auto average_size_estimate = unaryFun!fn(first_read);

        _appender.reserve(_size / average_size_estimate);
        _appender.put(first_read);

        while (total_size <= _size && !_range.empty) {
            auto read = _range.front;
            if (_split_by_ref && (read.ref_id != first_read.ref_id)) {
                break;
            }
            total_size += unaryFun!fn(read);
            _appender.put(read);
            _range.popFront();
        }
    }
}

/// Split range in chunks where total amount of memory consumed by all reads 
/// in the chunk is roughly chunk_size bytes.
///
/// Parameter $(D split_by_ref) specifies that each chunk should contain reads
/// aligned to the same reference. In most cases, this simplifies post-processing,
/// but in some cases this is not required, therefore it is optional.
auto chunksConsumingLessThan(R)(R reads, size_t size_in_bytes, bool split_by_ref=true) {
    return ReadRangeSplitter!(R, "a.size_in_bytes")(reads, size_in_bytes, split_by_ref);
}

/// Split range in chunks each containing no more than N reads
auto chunksOfSize(R)(R reads, size_t N, bool split_by_ref=true) {
    return ReadRangeSplitter!(R, "1")(reads, N, split_by_ref);
}
