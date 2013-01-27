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
