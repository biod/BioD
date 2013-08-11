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
module bio.core.bgzf.inputstream;

import bio.core.bgzf.block;
import bio.core.bgzf.blockrange;
import bio.core.bgzf.virtualoffset;
import bio.core.utils.file;

import std.stdio;
import std.stream : Stream, ReadException, WriteException, SeekException, SeekPos;
import std.range;
import std.array;
import std.traits;
import std.algorithm;
import std.parallelism;

/// Interface on top of Stream, which glues decompressed blocks
/// together and provides virtualTell() method for determining
/// virtual offset in the stream.
///
/// NOTE: the class has no virtualSeek method. The reason is
///       that an implementing class should have no direct access 
///       to any of underlying streams, (which provide chunks 
///       and their start virtual offsets), the reason being
///       the single responsibility principle. 
///       
///       If you need to do "virtualSeek", seek to coffset
///       in stream providing compressed blocks, create new
///       stream for decompressing, and skip uoffset bytes
///       in the decompressed stream. In general, you probably
///       will want to use different  decompression function
///       for each use case. (For serial iteration - without any 
///       caching, for random access - with caching; and maybe 
///       you'll want several types of caching.)
///
abstract class IChunkInputStream : Stream {
    /// Returns: current virtual offset
    ///
    /// (For details about what virtual offset is, 
    /// see bgzfrange module documentation.)
    ///
    abstract VirtualOffset virtualTell();

    /// Read a slice from chunk stream.
    abstract ubyte[] readSlice(size_t n);

    /// Average compression so far
    abstract float average_compression_ratio() @property const;

    /// Size of underlying BAM file (0 if not available).
    abstract size_t compressed_file_size() @property const;
}

/// Extends DecompressedBgzfBlock with skip_start and skip_end members.
struct AugmentedDecompressedBgzfBlock {
    DecompressedBgzfBlock block;
    alias block this;
    ushort skip_start; /// how many bytes to skip from start
    ushort skip_end;   /// how many bytes to skip from end
}

/// Convenience function for making decompressed blocks compatible with ChunkInputStream.
AugmentedDecompressedBgzfBlock makeAugmentedBlock(DecompressedBgzfBlock block) {
    return AugmentedDecompressedBgzfBlock(block, 0, 0);
}

/**
  Class for turning range of AugmentedDecompressedBgzfBlock objects
  into a proper InputStream

  NOTE: an assumption is made that calling empty() on range is cheap.
  However, front() gets called only once for each element. That fits
  the philosophy of std.algorithm.map. 
 */
final class ChunkInputStream(ChunkRange) : IChunkInputStream
{

    static assert(is(ElementType!ChunkRange == AugmentedDecompressedBgzfBlock));

    /// Construct class instance from range of chunks.
    this(ChunkRange range, ulong compressed_file_size=0) 
    {
        _range = range;
        readable = true;
        writeable = false;
        seekable = false;

        if (_range.empty) {
            setEOF();
        }

        _its_time_to_get_next_chunk = true; // defer getting first chunk

        _compressed_file_size = compressed_file_size;
    }

    /// Read a slice from chunk stream.
    ///
    /// Behaviour: when current chunk contains >= n bytes,
    ///            a chunk slice is returned. 
    ///            Otherwise, memory copying occurs.
    override ubyte[] readSlice(size_t n) {
        debug { stderr.writeln("[debug][ChunkInputStream] reading chunk of ", n, " bytes"); }
        if (_range.empty()) {
            _start_offset = _end_offset;
            _cur = 0;
            setEOF();
            return null;
        }

        if (_its_time_to_get_next_chunk) {
            setupStream();
        }

        if (_cur + n < _len) {
            // can return a slice, 
            // remaining in the current chunk
            auto slice = _buf[_cur .. _cur + n];
            _cur += n;
            return slice;
        } else if (_cur + n == _len) {
            // end of current chunk
            auto slice = _buf[_cur .. _len];
            _range.popFront();
            _its_time_to_get_next_chunk = true;
            return slice;
        } else {
            // this branch will be executed rarely
            // (wish the compiler could understand that...)
            auto slice = uninitializedArray!(ubyte[])(n);
            readExact(slice.ptr, n);
            return slice;
        }
    }

    override size_t readBlock(void* buffer, size_t size) {

        if (_range.empty()) {
            // no more chunks
            _start_offset = _end_offset;
            _cur = 0;
            setEOF();
            debug { stderr.writeln("[debug][ChunkInputStream] read ", 0, " bytes"); }
            return 0;
        }

        if (_its_time_to_get_next_chunk) {
            setupStream();
        }

        size_t read = readBlockHelper(buffer, size);

        if (read < size) {
            _range.popFront(); // get next chunk
            setupStream();
        }

        if (read == 0) { // try with next chunk
            read = readBlockHelper(buffer, size);
        }

        if (_cur == _len) { // end of current chunk
            _range.popFront();
            _its_time_to_get_next_chunk = true;
        }

        debug { stderr.writeln("[debug][ChunkInputStream] read ", read, " bytes"); }
        return read;
    }

    override size_t writeBlock(const void* buffer, size_t size) {
        throw new WriteException("Stream is not writeable");
    }

    override ulong seek(long offset, SeekPos whence) {
        throw new SeekException("Stream is not seekable");
    }

    /// Returns: current virtual offset
    override VirtualOffset virtualTell() {
        assert(_cur < (1<<16));
        if (_its_time_to_get_next_chunk) {
            setupStream();
        }
        return VirtualOffset(_start_offset, cast(ushort)_cur);
    }

    override float average_compression_ratio() @property const {
        if (_total_compressed == 0) return 0.0;
        return cast(float)_total_uncompressed/_total_compressed;
    }

    override ulong compressed_file_size() @property const {
        return _compressed_file_size;
    }

private:
    ChunkRange _range;
    typeof(_range.front.start_offset) _start_offset;
    typeof(_range.front.end_offset) _end_offset;
    typeof(_range.front.decompressed_data) _buf;

    ulong _compressed_file_size;

    ulong _len;  // current data length
    ulong _cur;  // current position

    ulong _total_compressed; // compressed bytes read so far
    ulong _total_uncompressed; // uncompressed size of blocks read so far

    bool _its_time_to_get_next_chunk;

    // Effect:
    //
    // _range is untouched (only empty() and front() are used)
    //
    // _buf is filled with the contents of _range.front.decompressed_data
    //
    // _cur, _start_offset_, end_offset, _total_uncompressed, and
    // _total_compressed variables are updated appropriately
    //
    // If _range.front.decompressed_data.length == 0, readEOF is set
    // (the reason is that this means we've just encountered a special
    // BGZF block, indicating end-of-file).
    //
    // _its_time_to_get_next_chunk is set back to false
    void setupStream() {

        _its_time_to_get_next_chunk = false;

        _cur = 0;

        // we don't rely on the caller that it will set _start_offset and
        // _end_offset appropriately
        if (_range.empty) {
            _start_offset = _end_offset;
            setEOF();
            return;
        }

        auto tmp = _range.front; /// _range might be lazy, 
                                 /// so extra front() calls
                                 /// can cost a lot

        debug {
            /*
            import std.stdio;
            if (tmp.skip_start != 0 || tmp.skip_end != 0) {
                writeln("reading partial decompressed block: ");
                writeln("   skip_start   = ", tmp.skip_start);
                writeln("   skip_end     = ", tmp.skip_end);
                writeln("   start_offset = ", tmp.start_offset);
                writeln("   end_offset   = ", tmp.end_offset);
            }
            */
        }

        _cur = tmp.skip_start;
        _start_offset = tmp.start_offset;
        _end_offset = tmp.end_offset;

        _total_compressed += _end_offset - _start_offset - tmp.skip_start - tmp.skip_end;

        _buf = tmp.decompressed_data[0 .. $ - tmp.skip_end];
        _len = _buf.length;

        _total_uncompressed += _len;

        if (_len == 0) {
            setEOF();
        }

        return;
    }

    version(development)
    {
        final void setEOF() {
            import std.stdio;
            std.stdio.stderr.writeln("[info][chunkinputstream] len == 0, start offset = ", _start_offset,
                                     ", end offset = ", _end_offset);
            readEOF = true;
        }
    } else {
        final void setEOF() { readEOF = true; }
    }
    
    size_t readBlockHelper(void* buffer, size_t size) {
        ubyte* cbuf = cast(ubyte*) buffer;
        if (size + _cur > _len)
          size = cast(size_t)(_len - _cur);
        ubyte[] ubuf = cast(ubyte[])_buf[cast(size_t)_cur .. cast(size_t)(_cur + size)];
        cbuf[0 .. size] = ubuf[];
        _cur += size;
        return size;
    }
}

/**
   Returns: input stream wrapping given range of memory chunks
 */
auto makeChunkInputStream(R)(R range, ulong compressed_file_size=0) 
    if (is(Unqual!(ElementType!R) == AugmentedDecompressedBgzfBlock))
{
    return new ChunkInputStream!R(range, compressed_file_size);
}

/// ditto
auto makeChunkInputStream(R)(R range, ulong compressed_file_size=0)
    if (is(Unqual!(ElementType!R) == DecompressedBgzfBlock))
{
    return makeChunkInputStream(map!makeAugmentedBlock(range), compressed_file_size);
}

/// Convenient decorator for input streams.
final class BgzfInputStream : IChunkInputStream {

    private IChunkInputStream _stream;

    override size_t readBlock(void* buffer, size_t size) {
        return _stream.readBlock(buffer, size);
    }

    override ulong seek(long offset, SeekPos whence) {
        throw new SeekException("Stream is not seekable");
    }

    override size_t writeBlock(const void* buf, size_t size) {
        throw new WriteException("Stream is not writeable");
    }

    /// Methods implementing $(D IChunkInputStream)
    override VirtualOffset virtualTell() {
        return _stream.virtualTell();
    }

    /// ditto
    override ubyte[] readSlice(size_t n) {
        return _stream.readSlice(n);
    }

    /// ditto
    override float average_compression_ratio() @property const {
        return _stream.average_compression_ratio;
    }

    /// ditto
    override size_t compressed_file_size() @property const {
        return _stream.compressed_file_size;
    }

    bool readEOF() @property {
        return _stream.readEOF;
    }

    /// Read BGZF blocks from supplied $(D stream) using $(D task_pool)
    /// to do parallel decompression
    this(File file, TaskPool task_pool=taskPool) {
        auto file_size = file.isSeekable ? file.size : 0;
        auto bgzf_range = BgzfRange(file);

        auto decompressed_blocks = task_pool.map!decompressBgzfBlock(bgzf_range, 24);
        _stream = makeChunkInputStream(decompressed_blocks, file_size);

        readable = true;
        writeable = false;
        seekable = false;
    }
}
