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
module bio.bam.bgzf.block;

import bio.bam.constants;

import std.array : uninitializedArray;
import std.conv;
import std.zlib : crc32, ZlibException;
import etc.c.zlib;
import std.exception;

/**
  Structure representing BGZF block.
  In general, users shouldn't use it, as it is EXTREMELY low-level.
 */
struct BgzfBlock {
    // field types are as in the SAM/BAM specification
    // ushort ~ uint16_t, char ~ uint8_t, uint ~ uint32_t

    public ulong start_offset; /// start offset in the file, in bytes

    /// end offset in the file, in bytes
    public ulong end_offset() @property const {
        return start_offset + bsize + 1;
    }

    public ushort bsize; /// total Block SIZE minus one

    public ushort cdata_size; /// compressed data size

    /// A buffer is used to reduce number of allocations.
    ///
    /// Its size is max(bsize + 1, input_size)
    /// Initially, it contains compressed data, but is rewritten
    /// during decompressBgzfBlock -- indeed, who cares about
    /// compressed data after it has been uncompressed?
    public ubyte[] _buffer = void;

    /// If block has been already decompressed, result is undefined.
    public ubyte[] compressed_data() @property {
        return _buffer[0 .. cast(size_t)cdata_size];
    }

    public uint crc32;
    public uint input_size; /// size of uncompressed data

    hash_t toHash() const pure @safe nothrow {
        // since the block can be either compressed or decompressed,
        // returning CRC sum is the easiest and safest thing to do
        return crc32;
    }

    bool opEquals(const ref BgzfBlock other) pure @safe nothrow {
        return crc32 == other.crc32;
    }

    int opCmp(const ref BgzfBlock other) const pure @safe nothrow {
        return crc32 < other.crc32 ? -1 :
               crc32 > other.crc32 ?  1 : 0;
    }
}

/**
  Struct representing decompressed BgzfBlock

  Start offset is needed to be able to tell current virtual offset,
  and yet be able to decompress blocks in parallel.
 */
struct DecompressedBgzfBlock {
    ulong start_offset;
    ulong end_offset;
    ubyte[] decompressed_data;
}

/// Function for BGZF block decompression.
/// Reuses buffer allocated for storing compressed data,
/// i.e. after execution buffer of the passed $(D block)
/// is overwritten with uncompressed data.
DecompressedBgzfBlock decompressBgzfBlock(BgzfBlock block) {

    if (block.input_size == 0) {
        return DecompressedBgzfBlock(block.start_offset, 
                                     block.start_offset + block.bsize + 1,
                                     cast(ubyte[])[]); // EOF marker
        // TODO: add check for correctness of EOF marker
    }

    int err = void;

    // allocate buffer on the stack
    ubyte[BGZF_MAX_BLOCK_SIZE] uncompressed_buf = void;

    // check that block follows BAM specification
    enforce(block.input_size <= BGZF_MAX_BLOCK_SIZE, 
            "Uncompressed block size must be within " ~ 
            to!string(BGZF_MAX_BLOCK_SIZE) ~ " bytes");

    // for convenience, provide a slice
    auto uncompressed = uncompressed_buf[0 .. block.input_size];

    // set input data
    etc.c.zlib.z_stream zs;
    zs.next_in = cast(typeof(zs.next_in))block.compressed_data;
    zs.avail_in = to!uint(block.compressed_data.length);

    err = etc.c.zlib.inflateInit2(&zs, /* winbits = */-15);
    if (err)
    {
        throw new ZlibException(err);
    }

    // uncompress it into a buffer on the stack
    zs.next_out = cast(typeof(zs.next_out))uncompressed_buf.ptr;
    zs.avail_out = block.input_size;

    err = etc.c.zlib.inflate(&zs, Z_FINISH);
    switch (err)
    {
        case Z_STREAM_END:
            assert(zs.total_out == block.input_size);
            err = etc.c.zlib.inflateEnd(&zs);
            if (err != Z_OK) {
                throw new ZlibException(err);
            }
            break;
        default:
            etc.c.zlib.inflateEnd(&zs);
            throw new ZlibException(err);
    }

    assert(block.crc32 == crc32(0, uncompressed[]));

    // Now copy back to block._buffer, overwriting existing data.
    // It should have enough bytes already allocated.
    assert(block._buffer.length >= block.input_size);
    block._buffer[0 .. block.input_size] = uncompressed[];

    return DecompressedBgzfBlock(block.start_offset, 
                                 block.start_offset + block.bsize + 1,
                                 block._buffer[0 .. block.input_size]);
}
