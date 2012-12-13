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
