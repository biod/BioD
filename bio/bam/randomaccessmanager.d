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

import bio.core.bgzf.blockrange;
import bio.core.bgzf.virtualoffset;
import bio.core.bgzf.inputstream;
import bio.core.utils.memoize;
import bio.core.utils.range;
import bio.core.utils.stream;

import std.system;
import std.algorithm;
import std.array;
import std.range;
import std.traits;
import std.exception;

import std.parallelism;

// keeps task pool together with block
struct BgzfBlockAux {
    TaskPool task_pool;
    BgzfBlock block;
    alias block this;

    hash_t toHash() const pure @safe nothrow { return block.toHash(); }

    bool opEquals(const ref BgzfBlockAux other) pure @safe nothrow {
        return block == other.block;
    }

    int opCmp(const ref BgzfBlockAux other) const pure @safe nothrow {
        return block.opCmp(other.block);
    }
}

// BgzfBlockAux -> Task
auto decompressTask(BgzfBlockAux b) {
    auto t = task!decompressBgzfBlock(b.block);
    b.task_pool.put(t);
    return t;
}

// BgzfBlockAux -> Task
private alias memoize!(decompressTask, 512, 
                       FifoCache, BgzfBlockAux) memDecompressTask;

// (BgzfBlock, TaskPool) -> DecompressedBgzfBlock
auto decompressSerial(BT)(BT block_and_pool) {
    return decompress(block_and_pool).yieldForce();
}

// (BgzfBlock, TaskPool) -> Task
auto decompress(BT)(BT block_and_pool) { 
    auto data = BgzfBlockAux(block_and_pool[1], block_and_pool[0]);
    return memDecompressTask(data);
}

// ([BgzfBlock], TaskPool) -> [DecompressedBgzfBlock]
auto parallelUnpack(BR)(BR bgzf_range, TaskPool pool, size_t n_threads = 0) {
    if (n_threads == 0)
        n_threads = max(pool.size(), 1);

    auto tasks = bgzf_range.zip(repeat(pool)).map!decompress();
    return tasks.prefetch(n_threads).map!"a.yieldForce()"();
}

debug {
    import std.stdio;
}

/// Class which random access tasks are delegated to.
class RandomAccessManager {

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

    /// Get new IChunkInputStream starting from specified virtual offset.
    IChunkInputStream createStreamStartingFrom(VirtualOffset offset, bool parallel=true) {

        auto _stream = new bio.core.utils.stream.File(_filename);
        auto _compressed_stream = new EndianStream(_stream, Endian.littleEndian);
        _compressed_stream.seekSet(cast(size_t)(offset.coffset));

        auto n_threads = parallel ? max(task_pool.size, 1) : 1;
        auto blocks = BgzfRange(_compressed_stream).parallelUnpack(task_pool, n_threads);

        static auto helper(R)(R decompressed_range, VirtualOffset offset) {

            auto adjusted_front = AugmentedDecompressedBgzfBlock(decompressed_range.front,
                                                                 offset.uoffset, 0); 
            decompressed_range.popFront();
            auto adjusted_range = chain(repeat(adjusted_front, 1), 
                                        map!makeAugmentedBlock(decompressed_range));

            return cast(IChunkInputStream)makeChunkInputStream(adjusted_range);
        }

        return helper(blocks, offset);
    }

    /// Get single read at a given virtual offset.
    /// Every time new stream is used.
    BamRead getReadAt(VirtualOffset offset) {
        auto stream = createStreamStartingFrom(offset);
        return bamReadRange(stream, _reader).front.dup;
    }

    /// Get BGZF block at a given offset.
    BgzfBlock getBgzfBlockAt(ulong offset) {
        auto stream = new bio.core.utils.stream.File(_filename);
        stream.seekSet(offset);
        return BgzfRange(stream).front;
    }

    /// Get reads between two virtual offsets. First virtual offset must point
    /// to a start of an alignment record.
    ///
    /// If $(D task_pool) is not null, it is used for parallel decompression. Otherwise, decompression is serial.
    auto getReadsBetween(VirtualOffset from, VirtualOffset to) {
        IChunkInputStream stream = createStreamStartingFrom(from);

        static bool offsetTooBig(BamReadBlock record, VirtualOffset vo) {
            return record.end_virtual_offset > vo;
        }

        return until!offsetTooBig(bamReadRange!withOffsets(stream, _reader), to);
    }

    bool found_index_file() @property {
        return _found_index_file;
    }
    private bool _found_index_file = false; // overwritten in constructor if filename is provided

    /// BAI file
    ref const(BaiFile) getBai() const {
        return _bai;
    }

    /// Get BAI chunks containing all alignment records overlapping specified region
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
    auto getReads(alias IteratePolicy=withOffsets)(int ref_id, uint beg, uint end) {
        auto chunks = array(nonOverlappingChunks(getChunks(ref_id, beg, end)));

        debug {
            /*
            import std.stdio;
            writeln("[random access] chunks:");
            writeln("    ", chunks);
            */
        }

        // General plan:
        //
        // chunk[0] -> bgzfRange[0] |
        // chunk[1] -> bgzfRange[1] | (2)
        //         ....             | -> (joiner(bgzfRange), [start/end v.o.])
        // chunk[k] -> bgzfRange[k] |                      |
        //         (1)                     /* parallel */  V                       (3)
        //                                  (unpacked blocks, [start/end v.o.])
        //                                                 |
        //                                                 V                       (4)
        //                                     (modified unpacked blocks)
        //                                                 |
        //                                                 V                       (5)
        //                                          IChunkInputStream
        //                                                 |
        //                                                 V                       (6)
        //                                 filter out non-overlapping records
        //                                                 |
        //                                                 V
        //                                              that's it!

        auto bgzf_range = getJoinedBgzfRange(chunks);                               // (2)
        auto decompressed_blocks = getUnpackedBlocks(bgzf_range, task_pool);        // (3)
        auto augmented_blocks = getAugmentedBlocks(decompressed_blocks, chunks);    // (4)
        IChunkInputStream stream = makeChunkInputStream(augmented_blocks);          // (5)
        auto reads = bamReadRange!IteratePolicy(stream, _reader);
        return filterBamReads(reads, ref_id, beg, end);                             // (6)
    }

private:
    
    string _filename;
    BaiFile _bai;
    BamReader _reader;
    TaskPool _task_pool;
    size_t _buffer_size;

    TaskPool task_pool() @property {
        if (_task_pool is null)
            _task_pool = taskPool;
        return _task_pool;
    }
        
public:

    // Let's implement the plan described above!

    // (1) : (Chunk, Stream) -> [BgzfBlock]
    static struct ChunkToBgzfRange {
        static bool offsetTooBig(BgzfBlock block, ulong offset) {
            return block.start_offset > offset;
        }

        private {
            Chunk _chunk;
            Stream _stream;
            bool _init = false;
            Until!(offsetTooBig, BgzfRange, ulong) _range;
        }

        this(Chunk chunk, Stream stream) {
            _chunk = chunk;
            _stream = stream;
        }

        auto front() @property { init(); return _range.front; }
        void popFront() { init(); _range.popFront(); }
        bool empty() @property { init(); return _range.empty; }

        private void init() {
            if (!_init) {
                _stream.seekSet(cast(size_t)_chunk.beg.coffset);
                _range = until!offsetTooBig(BgzfRange(_stream), 
                                            _chunk.end.coffset);
                _init = true;
            }
        }
    }

    // (2) : Chunk[] -> [BgzfBlock]
    auto getJoinedBgzfRange(Chunk[] bai_chunks) {
        Stream file = new bio.core.utils.stream.File(_filename);
        Stream stream = new BufferedStream(file, _buffer_size);

        ChunkToBgzfRange[] bgzf_ranges;
        bgzf_ranges.length = bai_chunks.length;
        foreach (i, ref range; bgzf_ranges) {
            range = ChunkToBgzfRange(bai_chunks[i], stream);
        }
        auto bgzf_blocks = joiner(bgzf_ranges);
        return bgzf_blocks;
    }

    // (3) : ([BgzfBlock], TaskPool) -> [DecompressedBgzfBlock]
    static auto getUnpackedBlocks(R)(R bgzf_range, TaskPool pool) {
        version(serial) {
            return bgzf_range.parallelUnpack(pool, 1);
        } else {
            return bgzf_range.parallelUnpack(pool);
        }
    }

    // (4) : ([DecompressedBgzfBlock], Chunk[]) -> [AugmentedDecompressedBgzfBlock]

    // decompressed blocks:
    // [.....][......][......][......][......][......][.....][....]
    //
    // what we need (chunks):
    //   [.........]  [.........]        [...........]  [..]
    //
    // Solution: augment decompressed blocks with skip_start and skip_end members
    //           and teach ChunkInputStream to deal with ranges of such blocks.
    static struct AugmentedBlockRange(R) {
        this(R blocks, Chunk[] bai_chunks) {
            _blocks = blocks;
            if (_blocks.empty) {
                _empty = true;
            } else {
                _cur_block = _blocks.front;
                _blocks.popFront();
            }
            _chunks = bai_chunks[];
        }

        bool empty() @property {
            return _empty;
        }

        AugmentedDecompressedBgzfBlock front() @property {
            AugmentedDecompressedBgzfBlock result;
            result.block = _cur_block;

            if (_chunks.empty) {
                return result;
            }

            if (beg.coffset == result.start_offset) {
                result.skip_start = beg.uoffset;
            }

            if (end.coffset == result.start_offset) {
                auto to_skip = result.decompressed_data.length - end.uoffset;
                assert(to_skip <= ushort.max);
                result.skip_end = cast(ushort)to_skip;
            }

            return result;
        }

        void popFront() {
            if (_cur_block.start_offset == end.coffset) {
                _chunks = _chunks[1 .. $];
            }
            if (_blocks.empty) {
                _empty = true;
                return;
            }
            _cur_block = _blocks.front;
            _blocks.popFront();
        }

        private {
            R _blocks;
            ElementType!R _cur_block;
            bool _empty;
            Chunk[] _chunks;

            VirtualOffset beg() @property {
                return _chunks[0].beg;
            }

            VirtualOffset end() @property {
                return _chunks[0].end;
            }
        }
    }

    static auto getAugmentedBlocks(R)(R decompressed_blocks, Chunk[] bai_chunks) {
        return AugmentedBlockRange!R(decompressed_blocks, bai_chunks);
    }

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
