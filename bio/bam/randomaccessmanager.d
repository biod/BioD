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
/**
  Module for random access operations on BAM file.
 */
module bio.bam.randomaccessmanager;

import bio.bam.constants;
import bio.bam.reader;
import bio.bam.read;
import bio.bam.readrange;
import bio.bam.baifile;
import bio.bam.bai.utils.algo;

import bio.core.bgzf.block;
import bio.core.bgzf.virtualoffset;
import bio.core.bgzf.inputstream;
import bio.core.bgzf.constants;
import bio.core.bgzf.chunk;
import bio.core.utils.range;
import bio.core.utils.stream;

import std.system;
import std.algorithm;
import std.array;
import std.range;
import std.traits;
import std.exception;

import std.parallelism;

debug {
    import std.stdio;
}

/// Class which random access tasks are delegated to.
class RandomAccessManager {

    void setCache(BgzfBlockCache cache) {
        _cache = cache;
    }
    
    void setTaskPool(TaskPool task_pool) {
        _task_pool = task_pool;
    }

    void setBufferSize(size_t buffer_size) {
        _buffer_size = buffer_size;
    }

    /// Constructs new manager for BAM file
    this(string filename) {
        _filename = filename;
    }

    /// ditto
    this(BamReader reader) {
        _reader = reader;
        _filename = reader.filename;
    }

    /// Constructs new manager with given index file.
    /// This allows to do random-access interval queries.
    ///
    /// Params:
    ///     filename =  location of BAM file
    ///     bai  =  index file
    this(string filename, ref BaiFile bai) {
        _filename = filename;
        _bai = bai;
        _found_index_file = true;
    }

    /// ditto
    this(BamReader reader, ref BaiFile bai) {
        _reader = reader;
        _filename = reader.filename;
        _bai = bai;
        _found_index_file = true;
    }

    /// If file ends with EOF block, returns virtual offset of the start of EOF block.
    /// Otherwise, returns virtual offset of the physical end of file.
    VirtualOffset eofVirtualOffset() const {
        ulong file_offset = std.file.getSize(_filename);
        if (hasEofBlock()) {
            return VirtualOffset(file_offset - BAM_EOF.length, 0);
        } else {
            return VirtualOffset(file_offset, 0);
        }
    }

    /// Returns true if the file ends with EOF block, and false otherwise.
    bool hasEofBlock() const {
        auto _stream = new bio.core.utils.stream.File(_filename);
        if (_stream.size < BAM_EOF.length) {
            return false;
        }

        ubyte[BAM_EOF.length] buf;
        _stream.seekEnd(-cast(int)BAM_EOF.length);

        _stream.readExact(&buf, BAM_EOF.length);
        if (buf != BAM_EOF) {
            return false;
        }

        return true;
    }

    /// Get new BgzfInputStream starting from specified virtual offset.
    BgzfInputStream createStreamStartingFrom(VirtualOffset offset)
    {
        auto _stream = new bio.core.utils.stream.File(_filename);
        auto _compressed_stream = new EndianStream(_stream, Endian.littleEndian);
        _compressed_stream.seekSet(cast(size_t)(offset.coffset));
        auto supplier = new StreamSupplier(_compressed_stream, offset.uoffset);
        auto bgzf_stream = new BgzfInputStream(supplier, _task_pool, _cache);
        return bgzf_stream;
    }

    /// Get single read at a given virtual offset.
    /// Every time new stream is used.
    BamRead getReadAt(VirtualOffset offset) {
        auto stream = createStreamStartingFrom(offset);
        
        bool old_mode = _reader._seqprocmode;
        _reader._seqprocmode = true;
        auto read = bamReadRange(stream, _reader).front.dup;
        _reader._seqprocmode = old_mode;
        return read;
    }

    /// Get BGZF block at a given offset.
    BgzfBlock getBgzfBlockAt(ulong offset) {
        auto fstream = new bio.core.utils.stream.File(_filename);
        auto stream = new EndianStream(fstream, Endian.littleEndian);
        stream.seekSet(offset);
        BgzfBlock block = void;
        ubyte[BGZF_MAX_BLOCK_SIZE] buf = void;
        fillBgzfBufferFromStream(stream, true, &block, buf.ptr);
        block._buffer = block._buffer.dup;
        return block;
    }

    /// Get reads between two virtual offsets. First virtual offset must point
    /// to a start of an alignment record.
    auto getReadsBetween(VirtualOffset from, VirtualOffset to) {
        auto stream = createStreamStartingFrom(from);

        static bool offsetTooBig(BamReadBlock record, VirtualOffset vo) {
            return record.end_virtual_offset > vo;
        }

        return until!offsetTooBig(bamReadRange!withOffsets(stream, _reader), to);
    }

    bool found_index_file() @property const {
        return _found_index_file;
    }
    private bool _found_index_file = false; // overwritten in constructor if filename is provided

    /// BAI file
    ref const(BaiFile) getBai() const {
        enforce(found_index_file, "BAM index file (.bai) must be provided");
        return _bai;
    }

    /// Get BAI chunks containing all alignment records overlapping the region
    Chunk[] getChunks(int ref_id, int beg, int end) {
        enforce(found_index_file, "BAM index file (.bai) must be provided");
        enforce(ref_id >= 0 && ref_id < _bai.indices.length, "Invalid reference sequence index");

        // Select all bins that overlap with [beg, end).
        // Then from such bins select all chunks that end to the right of min_offset.
        // Sort these chunks by leftmost coordinate and remove all overlaps.

        auto min_offset = _bai.indices[ref_id].getMinimumOffset(beg);

        Chunk[] bai_chunks;
        foreach (b; _bai.indices[ref_id].bins) {
            if (!b.canOverlapWith(beg, end))
                continue;

            foreach (chunk; b.chunks) {
                if (chunk.end > min_offset) {
                    bai_chunks ~= chunk;

                    // optimization
                    if (bai_chunks[$-1].beg < min_offset) {
                        bai_chunks[$-1].beg = min_offset;
                    }
                }
            }
        }

        sort(bai_chunks);

        return bai_chunks;
    }

    /// Fetch alignments with given reference sequence id, overlapping [beg..end)
    auto getReads(alias IteratePolicy=withOffsets)(int ref_id, uint beg, uint end)
    {
        auto chunks = array(nonOverlappingChunks(getChunks(ref_id, beg, end)));
        auto fstream = new bio.core.utils.stream.File(_filename);
        auto compressed_stream = new EndianStream(fstream, Endian.littleEndian);
        auto supplier = new StreamChunksSupplier(compressed_stream, chunks);
        auto stream = new BgzfInputStream(supplier, _task_pool, _cache);
        auto reads = bamReadRange!IteratePolicy(stream, _reader);
        return filterBamReads(reads, ref_id, beg, end);
    }

private:
    
    string _filename;
    BaiFile _bai;
    BamReader _reader;
    TaskPool _task_pool;
    size_t _buffer_size;

    BgzfBlockCache _cache;

    TaskPool task_pool() @property {
        if (_task_pool is null)
            _task_pool = taskPool;
        return _task_pool;
    }
        
public:

    static struct BamReadFilter(R) {
        this(R r, int ref_id, uint beg, uint end) {
            _range = r;
            _ref_id = ref_id;
            _beg = beg;
            _end = end;
            findNext();
        }

        bool empty() @property {
            return _empty;
        }

        ElementType!R front() @property {
            return _current_read;
        }
        
        void popFront() {
            _range.popFront();
            findNext();
        }

    private: 
        R _range;
        int _ref_id;
        uint _beg;
        uint _end;
        bool _empty;
        ElementType!R _current_read;

        void findNext() {
            if (_range.empty) {
                _empty = true;
                return;
            }
            while (!_range.empty) {
                _current_read = _range.front;

                // BamReads are sorted first by ref. ID.
                auto current_ref_id = _current_read.ref_id;
                if (current_ref_id > _ref_id) {
                    // no more records for this _ref_id
                    _empty = true;
                    return;
                } else if (current_ref_id < _ref_id) {
                    // skip reads referring to sequences
                    // with ID less than ours
                    _range.popFront();
                    continue;
                }

                if (_current_read.position >= _end) {
                    _empty = true;
                    // As reads are sorted by leftmost coordinate,
                    // all remaining alignments in _range 
                    // will not overlap the interval as well.
                    // 
                    //                  [-----)
                    //                  . [-----------)
                    //                  .  [---)
                    //                  .    [-------)
                    //                  .         [-)
                    //    [beg .....  end)
                    return;
                }

                if (_current_read.position > _beg) {
                    return; // definitely overlaps
                }

                if (_current_read.position +
                    _current_read.basesCovered() <= _beg) 
                {
                    /// ends before beginning of the region
                    ///  [-----------)
                    ///               [beg .......... end)
                    _range.popFront();
                    /// Zero-length reads are also considered non-overlapping,
                    /// so for consistency the inequality 12 lines above is strict.
                } else {
                    return; /// _current_read overlaps the region
                }
            }
            _empty = true; 
        }
    }

    // Get range of alignments sorted by leftmost coordinate,
    // together with an interval [beg, end),
    // and return another range of alignments which overlap the region.
    static auto filterBamReads(R)(R r, int ref_id, uint beg, uint end) 
    {
        return BamReadFilter!R(r, ref_id, beg, end);
    }
}
