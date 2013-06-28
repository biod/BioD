/*
    This file is part of BioD.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
/// Writing a script/tool for processing BAM data often starts this way:
/// 
/// ------------------------
/// import bio.bam.reader;
///
/// void main(string[] args) {
///     auto bam = new BamReader(args[1]); // open BAM file
///     foreach (read; bam.reads) {        // iterate through its reads
///         if (read.is_unmapped)
///             continue;                  // maybe skip unmapped ones
///         ...
///     }
/// }
/// ------------------------
/// 
/// Or, if a specific interval on the reference sequence is to be explored:
/// ------------------------
/// import bio.bam.pileup;
/// ...
/// auto reads = bam["chr7"][50_000 .. 60_000]; // BAI index is required
/// foreach (column; makePileup(reads)) { ... } // see $(PMODULE pileup) docs
/// ------------------------
module bio.bam.reader;

import bio.bam.abstractreader;
public import bio.sam.header;
public import bio.bam.reference;
public import bio.bam.read;
public import bio.bam.tagvalue;
public import bio.bam.readrange;
import bio.bam.randomaccessmanager;
import bio.bam.baifile;
import bio.bam.bai.indexing;
import bio.core.utils.range;
import bio.core.utils.stream;
public import bio.core.bgzf.blockrange;
public import bio.core.bgzf.inputstream;
public import bio.core.bgzf.virtualoffset;

import std.system;
import std.stdio;
import std.algorithm;
import std.range;
import std.conv;
import std.exception;
import std.parallelism;
import std.array;
import core.stdc.stdio;
import std.string;

/**
  BAM file reader, featuring parallel decompression of BGZF blocks.
 */
class BamReader : IBamSamReader {

    /**
      Creates reader associated with file or stream.
      (If stream constructor is used, no random access is possible.)
      $(BR)
      Optionally, task pool can be specified. 
      It will be used to unpack BGZF blocks in parallel.

      Example:
      -------------------------------------------
      import std.parallelism, bio.bam.reader;
      void main() {
        auto pool = new TaskPool(4); // use 4 threads
        scope (exit) pool.finish();  // don't forget!
        auto bam = new BamReader("file.bam", pool);
        ...
      }
      -------------------------------------------
     */
    this(std.stream.Stream stream, 
         std.parallelism.TaskPool task_pool = std.parallelism.taskPool) {
        _source_stream = new EndianStream(stream, Endian.littleEndian);
        _task_pool = task_pool;

        if (stream.seekable) {
            _stream_is_seekable = true;
        }

        initializeStreams();

        auto magic = _bam.readString(4);
        
        enforce(magic == "BAM\1", "Invalid file format: expected BAM\\1");

        readSamHeader();
        readReferenceSequencesInfo();

        // right after construction, we are at the beginning
        //                           of the list of reads

        if (_stream_is_seekable) {
            _reads_start_voffset = _decompressed_stream.virtualTell();
        }
    }

    /// ditto
    this(string filename, std.parallelism.TaskPool task_pool) {

        _filename = filename;
        _source_stream = getNativeEndianSourceStream();
        this(_source_stream, task_pool);
    }

    /// ditto
    this(string filename) {
        this(filename, std.parallelism.taskPool);
    }
 
    /**
      True if BAI file was found for this BAM file.
      This is necessary for any random-access operations.
      $(BR)
      Looks for files in the same directory which filename
      is either the file name of BAM file with '.bai' appended,
      or with the last extension replaced with '.bai'
      (that is, for $(I file.bam) paths $(I file.bai) and 
      $(I file.bam.bai) will be checked)
     */
    bool has_index() @property {
        return _random_access_manager.found_index_file;
    }

    /**
      Creates BAI file. If $(I overwrite) is false, it won't touch
      existing index if it is already found.
     */
    void createIndex(bool overwrite = false) {
        if (has_index && !overwrite)
            return;
        Stream stream = new BufferedFile(filename ~ ".bai", FileMode.OutNew);
        scope(exit) stream.close();
        bio.bam.bai.indexing.createIndex(this, stream);
        _bai_status = BaiStatus.notInitialized;
    }

    /** Filename, if the object was created via file name constructor,
        $(D null) otherwise.
     */
    string filename() @property const {
        return _filename;
    }

    /// If file ends with EOF block, returns virtual offset of the start of EOF block.
    /// Otherwise, returns virtual offset of the physical end of file.
    bio.core.bgzf.virtualoffset.VirtualOffset eofVirtualOffset() {
        return _random_access_manager.eofVirtualOffset();
    }

    /// Get BGZF block at a given file offset.
    bio.core.bgzf.block.BgzfBlock getBgzfBlockAt(ulong offset) {
        return _random_access_manager.getBgzfBlockAt(offset);
    }

    /**
      Returns: SAM header of the BAM file
     */
    bio.sam.header.SamHeader header() @property {
        if (_header is null) {
            synchronized {
                if (_header is null) {
                    _header = new SamHeader(_headertext);
                    _headertext = null;
                }
            }
        }
        return _header;
    }

    /**
        Returns: information about reference sequences
     */
    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const {
        return _reference_sequences;
    }

    /**
        Range of all alignment records in the file.
        $(BR)
        Element type of the returned range depends on the policy. 
        Default one is $(DPREF2 bam, readrange, withoutOffsets),
        in this case range element type is $(DPREF2 bam, read, BamRead).
        $(BR)
        The other option is $(DPREF2 bam, readrange, withOffsets),
        which allows to track read virtual offsets in the file.
        In this case range element type is $(DPREF2 bam, readrange, BamReadBlock).

        Example:
        ----------------------------------
        import bio.bam.readrange;
        ...
        auto bam = new BamReader("file.bam");
        auto reads = bam.reads!withOffsets();
        writeln(reads.front.start_virtual_offset);
        ----------------------------------
     */
    auto reads(alias IteratePolicy=bio.bam.readrange.withoutOffsets)() @property {
        auto _decompressed_stream = getDecompressedBamReadStream();
        return bamReadRange!IteratePolicy(_decompressed_stream, this);
    }

    static struct ReadsWithProgressResult(alias IteratePolicy, R, S) {
        this(R range, S stream, void delegate(lazy float p) progressBarFunc) {
            _range = range;
            _stream = stream;
            _progress_bar_func = progressBarFunc;
        }

        static if (__traits(identifier, IteratePolicy) == "withOffsets") {
            auto front() @property {
                return _range.front;
            } 
        } else static if (__traits(identifier, IteratePolicy) == "withoutOffsets") {
            auto front() @property {
                return _range.front.read;
            }
        } else static assert(0, __traits(identifier, IteratePolicy));

        bool empty() @property {
            return _range.empty;
        }

        void popFront() {
            _bytes_read += _range.front.read.size_in_bytes;
            _range.popFront();
            if (_progress_bar_func !is null) {
                _progress_bar_func(min(1.0, 
                    cast(float)_bytes_read / (_stream.compressed_file_size * 
                                              _stream.average_compression_ratio)));
            }
        }

        private R _range;
        private S _stream;
        private size_t _bytes_read;
        private void delegate(lazy float p) _progress_bar_func;
    }

    /**
        Returns: range of all reads in the file, calling $(I progressBarFunc)
                 for each read. 
        $(BR)
        $(I progressBarFunc) will be called
        each time next alignment is read, with the argument being a number from [0.0, 1.0],
        which is estimated progress percentage.
        $(BR)
        Notice that $(I progressBarFunc) takes $(D lazy) argument, 
        so that the number of relatively expensive float division operations
        can be controlled by user.

        Example:
        ------------------------------------
        import std.functional, std.stdio, bio.bam.reader;
        void progress(lazy float p) {
            static uint n;
            if (++n % 63 == 0) writeln(p); // prints progress after every 63 records
        }
        ...
        foreach (read; bam.readsWithProgress(toDelegate(&progress))) {
            ...
        }
        ------------------------------------
    */
    auto readsWithProgress(alias IteratePolicy=bio.bam.readrange.withoutOffsets)
        (void delegate(lazy float p) progressBarFunc) 
    {
        auto _decompressed_stream = getDecompressedBamReadStream();
        auto reads_with_offsets = bamReadRange!withOffsets(_decompressed_stream, this);
       
        alias ReadsWithProgressResult!(IteratePolicy, 
                       typeof(reads_with_offsets), IChunkInputStream) Result;
        
        return Result(reads_with_offsets, _decompressed_stream, progressBarFunc);
    }

    /// Part of IBamSamReader interface
    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads!withoutOffsets());
    }

    /**
      Returns: the read which starts at a given virtual offset.
     */
    bio.bam.read.BamRead getReadAt(bio.core.bgzf.virtualoffset.VirtualOffset offset) {
        enforce(_random_access_manager !is null);
        return _random_access_manager.getReadAt(offset);
    }

    /**
      Returns: all reads located between two virtual offsets in the BAM file.

      $(BR)
      First offset must point to the start of an alignment record,
      and be strictly less than the second one.
      $(BR)
      For decompression, the task pool specified at the construction is used.
     */ 
    auto getReadsBetween(bio.core.bgzf.virtualoffset.VirtualOffset from, 
                         bio.core.bgzf.virtualoffset.VirtualOffset to) {
        enforce(from < to, "First offset must be strictly less than second");
        enforce(_stream_is_seekable, "Stream is not seekable");
        
        return _random_access_manager.getReadsBetween(from, to);
    }

    /**
      Get BAI chunks containing all reads that overlap specified region.
     */
    bio.bam.bai.chunk.Chunk[] getChunks(int ref_id, int beg, int end) {
        enforce(_random_access_manager !is null);
        enforce(beg < end);

        return _random_access_manager.getChunks(ref_id, beg, end);
    }

    /**
      Returns reference sequence with id $(I ref_id).
     */
    bio.bam.reference.ReferenceSequence reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return ReferenceSequence(_random_access_manager, 
                                 ref_id,
                                 _reference_sequences[ref_id]);
    }

    /**
      Returns reference sequence named $(I ref_name).

      Example:
      ---------------------------
      import std.stdio, bio.bam.reader;
      ...
      auto bam = new BamReader("file.bam");
      writeln(bam["chr2"].length);
      ---------------------------
     */
    bio.bam.reference.ReferenceSequence opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " does not exist");
        auto ref_id = _reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    /**
      Check if reference named $(I ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        return null != (ref_name in _reference_sequence_dict);
    }

    /**
      Set buffer size for I/O operations. Values less than 4096 are disallowed.
      $(BR)
      This can help in multithreaded applications when several files are read 
      simultaneously (e.g. merging).
     */
    void setBufferSize(size_t buffer_size) {
        enforce(buffer_size >= 4096, "Buffer size must be >= 4096 bytes");
        _buffer_size = buffer_size;
        _random_access_manager.setBufferSize(buffer_size);
    }

private:
    
    string _filename;                       // filename (if available)
    Stream _source_stream;                  // compressed
    IChunkInputStream _decompressed_stream; // decompressed
    Stream _bam;                            // decompressed + endian conversion
    bool _stream_is_seekable;

    // Virtual offset at which alignment records start.
    VirtualOffset _reads_start_voffset;

    BaiFile _dont_access_me_directly_use_bai_file_for_that;
    enum BaiStatus {
        notInitialized,
        initialized,
        fileNotFound
    }
    BaiStatus _bai_status = BaiStatus.notInitialized;

    void initBai() {
        if (_bai_status == BaiStatus.notInitialized) {
            synchronized {
                try {
                    _dont_access_me_directly_use_bai_file_for_that = BaiFile(_filename);
                    _bai_status = BaiStatus.initialized;
                } catch (Exception e) {
                    _bai_status = BaiStatus.fileNotFound;
                }
            }
        }
    }

    // provides access to index file
    @property ref BaiFile _bai_file() { // initialized lazily
        initBai();
        return _dont_access_me_directly_use_bai_file_for_that;
    }; 

    RandomAccessManager _rndaccssmgr; // unreadable for a purpose
    @property RandomAccessManager _random_access_manager() {
        if (_rndaccssmgr is null) {
            synchronized {
                initBai();

                if (_bai_status == BaiStatus.initialized) {
                    _rndaccssmgr = new RandomAccessManager(this, _bai_file);
                } else {
                    _rndaccssmgr = new RandomAccessManager(this);
                }

                _rndaccssmgr.setTaskPool(_task_pool);
                _rndaccssmgr.setBufferSize(_buffer_size);
            }
        }
        return _rndaccssmgr;
    }

    SamHeader _header;
    string _headertext; // for lazy SAM header parsing
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict; /// name -> index mapping

    TaskPool _task_pool;
    size_t _buffer_size = 8192; // buffer size to be used for I/O

    Stream getNativeEndianSourceStream() {
        assert(_filename !is null);
        Stream file = new bio.core.utils.stream.File(_filename);
        return new BufferedStream(file, _buffer_size);
    }

    Stream getSeekableCompressedStream() {
        if (_stream_is_seekable) {
            if (_filename !is null) {
                auto file = getNativeEndianSourceStream();
                version(development)
                {
                    std.stdio.stderr.writeln("[info] file size: ", file.size);
                }
                return new EndianStream(file, Endian.littleEndian);
            } else {
                _source_stream.seekSet(0);
                return _source_stream;
            } 
        } else {
            return null;
        }
    }

    // get decompressed stream out of compressed BAM file
    IChunkInputStream getDecompressedStream() {

        auto compressed_stream = getSeekableCompressedStream();

        auto bgzf_range = (compressed_stream is null) ? BgzfRange(_source_stream) :
                                                        BgzfRange(compressed_stream);
        version(serial) {
            auto chunk_range = map!decompressBgzfBlock(bgzf_range);
        } else {
            auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 24);
        }

        if (compressed_stream !is null) {
            return makeChunkInputStream(chunk_range, cast(size_t)compressed_stream.size);
        } else {
            return makeChunkInputStream(chunk_range);
        }
    }


    // get decompressed stream starting from the first alignment record
    IChunkInputStream getDecompressedBamReadStream() {
        auto compressed_stream = getSeekableCompressedStream();

        if (compressed_stream !is null) {
            enforce(_reads_start_voffset != 0UL);

            compressed_stream.seekCur(_reads_start_voffset.coffset);
            auto bgzf_range = BgzfRange(compressed_stream);

            version(serial) {
                auto chunk_range = map!decompressBgzfBlock(bgzf_range);
            } else {
                auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 24);
            }
        
            auto sz = compressed_stream.size;
            auto stream = makeChunkInputStream(chunk_range, cast(size_t)sz);
            stream.readString(_reads_start_voffset.uoffset);
            return stream;
        } else {
            // must be initialized in initializeStreams()
            return _decompressed_stream;
        }
    }

    // sets up the streams and ranges
    void initializeStreams() {
        
        _decompressed_stream = getDecompressedStream();
        _bam = new EndianStream(_decompressed_stream, Endian.littleEndian); 
    }

    // initializes _header
    void readSamHeader() {
        int l_text;
        _bam.read(l_text);

        _headertext = to!string(_bam.readString(l_text));
    }

    // initialize _reference_sequences
    void readReferenceSequencesInfo() {
        int n_ref;
        _bam.read(n_ref);
        _reference_sequences = new ReferenceSequenceInfo[n_ref];
        foreach (i; 0..n_ref) {
            _reference_sequences[i] = ReferenceSequenceInfo(_bam);

            // provide mapping Name -> Index
            _reference_sequence_dict[_reference_sequences[i].name] = i;
        }
    }
}
