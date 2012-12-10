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
module bio.bam.readrange;

import bio.bam.read;
import bio.core.bgzf.inputstream;
import bio.core.bgzf.virtualoffset;
import bio.core.utils.switchendianness;

import std.stream;
import std.algorithm;
import std.system;

/// Tuple of virtual offset of the read, and the read itself.
struct BamReadBlock {
    VirtualOffset start_virtual_offset;
    VirtualOffset end_virtual_offset;
    BamRead read;
    alias read this;

    BamReadBlock dup() @property const {
        return BamReadBlock(start_virtual_offset, end_virtual_offset, read.dup);
    }
}

/// Policies for bamReadRange
mixin template withOffsets() {
    /**
        Returns: virtual offsets of beginning and end of the current read
                 plus the current read itself.
     */
    BamReadBlock front() @property {
        return BamReadBlock(_start_voffset, 
                            _stream.virtualTell(),
                            _current_record);
    }

    private VirtualOffset _start_voffset;

    private void beforeNextBamReadLoad() {
        _start_voffset = _stream.virtualTell();
    }
}

/// ditto
mixin template withoutOffsets() {
    /**
        Returns: current bamRead
     */
    ref BamRead front() @property {
        return _current_record;
    }

    private void beforeNextBamReadLoad() {}
}

class BamReadRange(alias IteratePolicy) 
{ 

    /// Create new range from BgzfInputStream.
    this(ref BgzfInputStream stream) {
        _stream = stream;
        _endian_stream = new EndianStream(_stream, Endian.littleEndian);
        readNext();
    }

    bool empty() @property const {
        return _empty;
    }

    mixin IteratePolicy;
    
    void popFront() {
        readNext();
    }

private:
    BgzfInputStream _stream;
    EndianStream _endian_stream;

    BamRead _current_record;
    bool _empty = false;

    /**
      Reads next bamRead block from stream.
     */
    void readNext() {

        // In fact, on BAM files containing a special EOF BGZF block
        // this condition will be always false!
        //
        // The reason is that we don't want to unpack next block just
        // in order to see if it's an EOF one or not.
        if (_stream.eof()) {
            _empty = true;
            return;
        }
     
        // In order to get the right virtual offset, we need to do it here.
        beforeNextBamReadLoad();

        // Here's where _empty is really set!
        int block_size = void;
        ubyte* ptr = cast(ubyte*)(&block_size);
        auto _read = 0;
        while (_read < int.sizeof) {
            auto _actually_read = _endian_stream.readBlock(ptr, int.sizeof - _read);
            if (_actually_read == 0) {
                version(development) {
                    import std.stdio;
                    stderr.writeln("[info][bamRead range] empty, read ", _read, " bytes, expected ", int.sizeof);
                }
                _empty = true;
                return;
            }
            _read += _actually_read;
            ptr += _actually_read;
        } 

        if (std.system.endian != Endian.littleEndian) {
            switchEndianness(&block_size, int.sizeof);
        }

        _current_record = BamRead(_stream.readSlice(block_size));
    }
}

/// Returns: lazy range of BamRead structs constructed from a given stream.
auto bamReadRange(alias IteratePolicy=withoutOffsets)(ref BgzfInputStream stream) {
    return new BamReadRange!IteratePolicy(stream);
}
