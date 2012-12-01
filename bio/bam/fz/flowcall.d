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
module bio.bam.fz.flowcall;

import bio.bam.tagvalue;
import bio.bam.fz.flowindex;

import bio.core.base;
import bio.core.utils.range;

import std.array;
import std.typecons;
import std.range;
import std.algorithm;

/// Flow base call
struct FlowCall {
    private {
        ushort _signal_intensity;
        Base _base;
    }

    /// Nucleotide
    Base base() @property const {
        return _base;
    }

    /// Signal intensity, normalized to homopolymer lengths
    float intensity() @property const {
        return _signal_intensity / 100.0;
    }

    /// round(intensity * 100.0)
    /// More efficient, because this is how intensities are stored in FZ tag.
    ushort intensity_value() @property const {
        return _signal_intensity;
    }
}

/// Flow call associated with a read
struct ReadFlowCall {
    private {
        ushort _signal_intensity;
        ushort _offset;
        ushort _called_len;
        Base _base;
        size_t _flow_index;
    }

    /// Called nucleotide
    Base base() @property const {
        return _base;
    }

    /// Set base to its complement
    void complement() {
        _base = _base.complement;
    }

    /// Called homopolymer length
    ushort length() @property const {
        return _called_len;
    }

    /// Zero-based position of the first nucleotide in the run,      
    /// relative to start of the read. Takes strandness into account.
    ushort offset() @property const {
        return _offset;
    }

    /// Signal intensity, normalized to homopolymer lengths
    float intensity() @property const {
        return _signal_intensity / 100.0;
    }

    /// round(intensity * 100.0)
    /// More efficient, because this is how intensities are stored in FZ tag.
    ushort intensity_value() @property const {
        return _signal_intensity;
    }

    /// Flow index (0-based)
    size_t flow_index() @property const {
        return _flow_index;
    }
}

/// Get flow calls from signal intensities and flow order.
auto flowCalls(ushort[] intensities, string flow_order) {
    
    static FlowCall flowCall(T)(T call) {
        return FlowCall(call[0], Base(call[1]));
    }

    return map!flowCall(zip(intensities, flow_order));
}

struct ReadFlowCallRange(S) 
    if (!is(S == class))
{
    private {
        string _flow_order = void;
        ushort[] _intensities = void;
        S _sequence = void;

        int _zf = void;
        Base _current_base = void;
        ushort _current_length = void;
        size_t _current_flow_index;
        ushort _current_offset;

        ushort _overlap;

        bool _empty;

        // consumes next homopolymer from the sequence,
        // and updates _current_base, _current_flow_index, 
        // _current_length appropriately
        void _doSetup() {
            if (_sequence.empty) {
                _empty = true;
                return;
            }

            // setup current base and current length
            _current_base = _sequence.front;
            _sequence.popFront();
            _current_length = 1;
            while (!_sequence.empty && _sequence.front == _current_base) {
                _sequence.popFront();
                ++_current_length;
            }

            // setup current flow index
            for ( ; _current_flow_index < _flow_order.length; ++_current_flow_index) {
                if (_flow_order[_current_flow_index] == _current_base) {
                    break;
                }
            }
        }
    }

    this(S seq, ushort[] intensities, string flow_order, ushort first_base_overlap, int zf) {
        _sequence = seq;
        _intensities = intensities;
        _flow_order = flow_order;
        _zf = zf;
        _overlap = first_base_overlap;

        if (_sequence.empty) {
            _empty = true;
        } else {
            _doSetup();
        }
    }

    bool empty() @property const {
        return _empty;
    }

    ReadFlowCall front() @property const {
        auto intensity = cast(ushort)(_intensities[_current_flow_index] - _overlap);
        return ReadFlowCall(intensity, _current_offset, _current_length, 
                            _current_base, _current_flow_index + _zf);
    }

    void popFront() {
        _current_offset += _current_length;

        ++_current_flow_index;
        _overlap = 0; // after first base it is always zero

        _doSetup();
    }

    ReadFlowCallRange!S save() @property {
        // bitwise copy
        // FIXME: is it safe?
        ReadFlowCallRange!S r = this;
        return r;
    }
}

private ReadFlowCallRange!S readFlowCallRange(S)(S seq, ushort[] intensities, 
                                                 string flow_order, ushort overlap, int zf)
{
    return ReadFlowCallRange!S(seq, intensities, flow_order, overlap, zf);
}


/// Get read flow calls. Takes ZF tag and strandness into account.
///
/// Tag name is an optional argument because it is not standard and will likely
/// be changed in the future (there was a proposal on samtools mailing list
/// to introduce standard FB tag).
ForwardRange!ReadFlowCall readFlowCalls(R)(R read, string flow_order, string key_sequence, string tag="ZF") {

    auto zf = cast(int)read[tag];
    Value fz_value = read["FZ"];
    ushort[] fz = *cast(ushort[]*)(&fz_value);

    flow_order = flow_order[zf .. $];
    auto intensities = fz[zf .. $];

    // key sequence is required because its last base can overlap with first called base
    ushort overlap = 0;

    Base5 base = read.is_reverse_strand ? read.sequence.back.complement : read.sequence.front;
    foreach_reverse (c; key_sequence) {
        if (c != base)
            break;
        overlap += 100;
    }

    if (!read.is_reverse_strand) {
        auto seq = read.sequence;
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order, overlap, zf));
    } else {
        auto seq = retro(map!"a.complement"(read.sequence));
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order, overlap, zf));
    }
}
