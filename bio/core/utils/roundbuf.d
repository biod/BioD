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
module bio.core.utils.roundbuf;

import std.exception;

/// Cyclic buffer
struct RoundBuf(T) {
   
    private {
        T[] _items = void;
        size_t _put;
        size_t _taken;
    }

    /** initializes round buffer of size $(D n) */
    this(size_t n) {
        _items = new T[n];
    }

    /// Input range primitives
    bool empty() @property const {
        return _put == _taken;
    }

    /// ditto
    auto ref front() @property {
        enforce(!empty, "buffer is empty");
        return _items[_taken % $];
    }

    /// ditto
    void popFront() {
        ++_taken;
    }

    /// Output range primitive
    void put(T item) {
        enforce(!full, "buffer is full");
        _items[_put % $] = item;
        ++_put;
    }

    /// Check if buffer is full
    bool full() @property const {
        return _put == _taken + _items.length;
    }
   
    /// Current number of elements
    size_t length() @property const {
        return _put - _taken;
    }
}

unittest {
    auto buf = RoundBuf!int(4);
    assert(buf.empty);

    buf.put(1);
    buf.put(2);
    assert(buf.length == 2);
    assert(buf.front == 1);
    buf.popFront();
    buf.put(1);
    buf.put(0);
    buf.put(3);
    assert(buf.full);
    buf.popFront();
    buf.put(4);
    buf.popFront();
    buf.popFront();
    assert(buf.front == 3);
    buf.popFront();
    assert(buf.front == 4);
}
