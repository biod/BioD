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
module bio.bam.utils.tagstoragebuilder;

import bio.bam.tagvalue;
import bio.bam.utils.value;

import std.array;

/// Struct for building tag storage, effectively managing memory.
struct TagStorageBuilder {
    private Appender!(ubyte[]) _storage;

    /// Return tag data (little-endian, in BAM format)
    ubyte[] data() @property {
        return _storage.data();
    }

    static TagStorageBuilder create() {
        TagStorageBuilder builder;
        builder._storage = appender!(ubyte[])();
        return builder;
    }

    /// Clear underlying storage
    void clear() {
        _storage.clear();
    }

    private void putByRef(string tag, ref Value v) {
        _storage.put(cast(ubyte[])tag);
        emplaceValue(_storage, v);
    }

    /// Append tag value to the storage
    void put(string tag, ref Value value) {
        putByRef(tag, value);
    }

    /// ditto
    void put(string tag, Value value) {
        putByRef(tag, value);
    }
}
