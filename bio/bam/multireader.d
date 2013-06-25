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
module bio.bam.multireader;

import bio.sam.header;
import bio.bam.reader;
import bio.bam.read;
import bio.bam.referenceinfo;
import bio.bam.utils.samheadermerger;

import std.algorithm;
import std.range;
import std.conv;
import std.parallelism;
import std.array;
import std.exception;
import std.typecons;

alias size_t FileId;

// ((BamRead, SamHeaderMerger), (BamRead, SamHeaderMerger)) -> bool
bool compare(T)(const auto ref T r1, const auto ref T r2) {
    assert(r1[1] == r2[1]);
    auto sorting_order = r1[1].merged_header.sorting_order;
    if (sorting_order == SortingOrder.coordinate)
        return compareCoordinates(r1[0], r2[0]);
    else if (sorting_order == SortingOrder.queryname)
        return compareReadNames(r1[0], r2[0]);
    else
        assert(0);
}

// (BamReader, SamHeaderMerger, FileId) -> [(BamRead, SamHeaderMerger, FileId)]
auto readRange(BamReader reader, SamHeaderMerger merger, FileId index) {
    return zip(reader.reads, repeat(merger), repeat(index));
}

// (BamReader, SamHeaderMerger, FileId, int, uint, uint) -> 
//                                    [(BamRead, SamHeaderMerger, FileId)]
auto readRange(BamReader reader, SamHeaderMerger merger, FileId index,
               int ref_id, uint start, uint end) 
{
    auto old_ref_id = cast(int)merger.ref_id_reverse_map[index][ref_id];
    auto reads = reader.reference(old_ref_id)[start .. end];
    return zip(reads, repeat(merger), repeat(index));
}

// ([BamReader], SamHeaderMerger) -> [[(BamRead, SamHeaderMerger, FileId)]]
auto readRanges(BamReader[] readers, SamHeaderMerger merger) {
    return readers.zip(repeat(merger), iota(readers.length))
                  .map!(t => readRange(t[0], t[1], t[2]))();
}

// ([BamReader], SamHeaderMerger, int, uint, uint) -> 
//                                    [[BamRead, SamHeaderMerger, FileId)]]
auto readRanges(BamReader[] readers, SamHeaderMerger merger, 
                int ref_id, uint start, uint end) 
{
    return readers.zip(repeat(merger), iota(readers.length), 
                       repeat(ref_id), repeat(start), repeat(end))
                  .map!(t => readRange(t[0], t[1], t[2], t[3], t[4], t[5]))();
}

// tweaks RG and PG tags, and reference sequence ID
// [[(BamRead, SamHeaderMerger, size_t)]] -> [[BamRead]]
auto adjustTags(R)(R reads_with_aux_info, TaskPool pool, size_t bufsize) 
    if (isInputRange!R) 
{
    static auto adj(U)(U tpl) { return adjustTags1(tpl[0], tpl[1], tpl[2]); }
    return reads_with_aux_info.zip(repeat(pool), repeat(bufsize))
                              .map!adj().array();
}

// [(BamRead, SamHeaderMerger, size_t)] -> [BamRead]
auto adjustTags1(R)(R reads_with_aux_info, TaskPool pool, size_t bufsize) 
    if (isInputRange!R) 
{
    return pool.map!adjustTags(reads_with_aux_info, bufsize);
}

// (BamRead, SamHeaderMerger, size_t) -> (BamRead, SamHeaderMerger)
auto adjustTags(R)(R read_with_aux_info) if (!isInputRange!R) {
    auto read = read_with_aux_info[0];
    auto merger = read_with_aux_info[1];
    auto file_id = read_with_aux_info[2];

    with (merger) {
        assert(file_id < ref_id_map.length);

        auto old_ref_id = read.ref_id;
        if (old_ref_id != -1 && old_ref_id in ref_id_map[file_id]) {
            auto new_ref_id = to!int(ref_id_map[file_id][old_ref_id]);
            if (new_ref_id != old_ref_id)
                read.ref_id = new_ref_id;
        } 

        auto program = read["PG"];
        if (!program.is_nothing) {
            auto pg_str = *(cast(string*)(&program));
            if (pg_str in program_id_map[file_id]) {
                auto new_pg = program_id_map[file_id][pg_str];
                if (new_pg != pg_str)
                    read["PG"] = new_pg;
            }
        }

        auto read_group = read["RG"];
        if (!read_group.is_nothing) {
            auto rg_str = *(cast(string*)(&read_group));
            if (rg_str in readgroup_id_map[file_id]) {
                auto new_rg = readgroup_id_map[file_id][rg_str];
                if (new_rg != rg_str)
                    read["RG"] = new_rg;
            }
        }
    }
    return tuple(read, merger);
}

///
class MultiBamReader {
  
    ///
    this(BamReader[] readers) {
        _readers = readers;
        _merger = new SamHeaderMerger(readers.map!(r => r.header)().array());

        auto n_references = _merger.merged_header.sequences.length;
        _reference_sequences = new ReferenceSequenceInfo[n_references];
        size_t i;
        foreach (line; _merger.merged_header.sequences) {
            _reference_sequences[i] = ReferenceSequenceInfo(line.name, line.length);
            _reference_sequence_dict[line.name] = i++;
        } 

        // TODO: maybe try to guess optimal size, based on the number of files?
        setBufferSize(1_048_576);
    }

    ///
    this(string[] filenames) {
        this(filenames.map!(fn => new BamReader(fn))().array());
    }

    ///
    this(string[] filenames, std.parallelism.TaskPool task_pool = taskPool) {
        this(filenames.zip(repeat(task_pool))
                      .map!(fn => new BamReader(fn[0], fn[1]))().array());
    }

    ///
    BamReader[] readers() @property {
        return _readers;
    }

    ///
    SamHeader header() @property {
        return _merger.merged_header;
    }

    ///
    auto reads() @property {
        return readRanges(_readers, _merger).adjustTags(task_pool, _adj_bufsz)
                                            .nWayUnion!compare().map!"a[0]"();
    }

    ///
    const(ReferenceSequenceInfo)[] reference_sequences() @property const {
        return _reference_sequences;
    }

    /**
      Check if reference named $(I ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        return null != (ref_name in _reference_sequence_dict);
    }

    /**
      Returns reference sequence with id $(I ref_id).
     */
    MultiBamReference reference(int ref_id) {
        enforce(ref_id >= 0, "Invalid reference index");
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return MultiBamReference(_readers, _merger, task_pool, _adj_bufsz,
                                 _reference_sequences[ref_id], ref_id);
    }

    /**
      Returns reference sequence named $(I ref_name).
     */
    MultiBamReference opIndex(string ref_name) {
        enforce(hasReference(ref_name), 
                "Reference with name " ~ ref_name ~ " does not exist");
        auto ref_id = cast(int)_reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    /// Sets buffer size for all readers (default is 1MB)
    void setBufferSize(size_t bytes) {
        foreach (reader; _readers)
            reader.setBufferSize(bytes);
    }

    private {
        BamReader[] _readers;
        SamHeaderMerger _merger;
        ReferenceSequenceInfo[] _reference_sequences;
        size_t[string] _reference_sequence_dict;
        TaskPool _task_pool;
        TaskPool task_pool() @property {
            if (_task_pool is null)
                _task_pool = taskPool;
            return _task_pool;
        }

        size_t _adj_bufsz = 512;
    }
}

///
struct MultiBamReference {
    private {
        BamReader[] _readers;
        SamHeaderMerger _merger;
        int _ref_id;
        ReferenceSequenceInfo _info;
        TaskPool _pool;
        size_t _adj_bufsz;
    }

    this(BamReader[] readers, SamHeaderMerger merger, 
         TaskPool task_pool, size_t adj_bufsize,
         ReferenceSequenceInfo info, int ref_id) 
    {
        _readers = readers;
        _merger = merger;
        _pool = task_pool;
        _adj_bufsz = adj_bufsize;
        _ref_id = ref_id;
        _info = info;
    }

    ///
    string name() @property const { return _info.name; }

    ///
    int length() @property const { return _info.length; }

    ///
    int ref_id() @property const { return _ref_id; }

    /// Get alignments overlapping [start, end) region.
    /// $(BR)
    /// Coordinates are 0-based.
    auto opSlice(uint start, uint end) {
        enforce(start < end, "start must be less than end");
        auto ranges = readRanges(_readers, _merger, ref_id, start, end);
        return ranges.adjustTags(_pool, _adj_bufsz)
                     .nWayUnion!compare().map!"a[0]"();
    }

    ///
    auto opSlice() {
        return opSlice(0, length);
    }
}
