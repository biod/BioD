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
module bio.bam.bai.utils.algo;

import bio.core.bgzf.chunk;

import std.range;
import std.algorithm;
import std.array;

struct NonOverlappingChunks(R) {

    this(R r) {
        _range = r;
        popFront();
    }

    bool empty() @property {
        return _empty;
    }

    auto front() @property {
        return _front;
    }

    void popFront() {
        if (!_range.empty()) {
            _front = _range.front;
            tryToJoinWithNextChunks();
        } else {
            _empty = true;
        }
    }

private:

    R _range;

    void tryToJoinWithNextChunks() {
        _range.popFront();
        while (!_range.empty()) {
            /// Get next element
            auto next = _range.front;
            /// It's presumed that chunks are sorted
            assert(next >= _front);
            /// Check if the chunk can be joined with the previous one
            if (_front.end >= next.beg) {
                /// update end of _front
                _front.end = max(_front.end, next.end);
                _range.popFront(); /// next is consumed
            } else {
                /// can't join
                break;
            }
        }
    }

    Chunk _front;
    bool _empty = false;
}

/// Params: 
///      r - range of chunks, sorted by leftmost coordinate
/// Returns:
///      range of non-overlapping chunks, covering the same subset
///      as original chunks
NonOverlappingChunks!R nonOverlappingChunks(R)(R r) 
    if (is(ElementType!R == Chunk)) 
{
    return NonOverlappingChunks!R(r);
}
