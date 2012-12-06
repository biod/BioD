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
module bio.bam.baseinfo;

import bio.core.base;

import bio.bam.read;
import bio.bam.fz.flowcall;

import std.range;
import std.conv;
import std.traits;
import std.typetuple;

///
struct MixinArg(T, string Tag) {
    T value;
    alias value this;
    alias Tag TagName;
}

/// Wrapper for arguments to $(D basesWith) function (see below).
/// Required to distinguish to which tag each parameter refers.
MixinArg!(T, Tag) arg(string Tag, T)(T value) {
    return MixinArg!(T, Tag)(value);
}

struct PerBaseInfo(R, Tags...) {

    private alias TypeTuple!("CIGAR", Tags) Extensions;

    // /////////////////////////////////////////////////////////////////////////
    //
    // Each 'extension' is a template with name TAGbaseInfo, containing 
    // a couple of mixin templates:
    // 
    // * resultProperties
    //      These are additional properties provided by the template
    //
    // * rangeMethods
    //      These describe how to proceed to the next base.
    //      The following methods must be implemented:
    //      
    //      - void setup(Args...)(const ref R read, Args args);
    //          Gets called during range construction. All constructor
    //          arguments are forwarded, and it's this function which
    //          is responsible for getting required parameters for this
    //          particular template.
    //
    //      - void populate(Result)(ref Result r);
    //          Populates fields of the result declared in resultProperties.
    //          Should run in O(1), just copying a few variables.
    //
    //      - void update(const ref R read);
    //          Encapsulates logic of moving to the next base and updating
    //          mixin variables correspondingly.
    //
    //      - void copy(Range)(const ref Range source, ref Range target);
    //          Gets called during $(D source.save). Therefore, any ranges
    //          used in mixin templates must be saved as well at that time.
    //
    // /////////////////////////////////////////////////////////////////////////

    private static string getResultProperties(Exts...)() {
        char[] result;
        foreach (ext; Exts) 
            result ~= "mixin " ~ ext ~ "baseInfo!R.resultProperties;".dup;
        return cast(string)result;
    }

    static struct Result {
        Base base;
        alias base this;

        string opCast(T)() if (is(T == string))
        {
            return to!string(base);
        }

        mixin(getResultProperties!Extensions());
    }

    private static string getRangeMethods(Exts...)() {
        char[] result;
        foreach (ext; Exts)
            result ~= "mixin " ~ ext ~ "baseInfo!R.rangeMethods " ~ ext ~ ";".dup;
        return cast(string)result;
    }

    mixin(getRangeMethods!Extensions());

    private void setup(string tag, Args...)(R read, Args args) {
        mixin(tag ~ ".setup(read, args);");
    }

    private void populate(string tag)(ref Result r) {
        mixin(tag ~ ".populate(r);");
    }

    private void update(string tag)() {
        mixin(tag ~ ".update(_read);");
    }

    private void copy(string tag)(ref typeof(this) other) {
        mixin(tag ~ ".copy(this, other);");
    }

    this(Args...)(R read, Args args) {
        _read = read;
        _seq = read.sequence;
        _rev = read.is_reverse_strand;

        foreach (t; Extensions) {
            setup!t(read, args);
        }
    }

    bool empty() @property const {
        return _seq.empty;
    }

    Result front() @property {
        Result r = void;
        r.base = _rev ? _seq.back.complement : _seq.front;
        foreach (t; Extensions)
            populate!t(r);
        return r;
    }

    void popFront() {
        moveToNextBase();
    }

    PerBaseInfo!(R, Tags) save() @property {
        PerBaseInfo!(R, Tags) r = void;
        r._read = _read.dup;
        r._seq = r._read.sequence;
        r._rev = _rev;
        foreach (t; Extensions)
            copy!t(r);
        return r;
    }

    private void moveToNextBase() {

        foreach (t; Extensions) {
            update!t();
        }

        if (_rev) 
            _seq.popBack();
        else 
            _seq.popFront();
    }

    private {
        R _read = void;
        typeof(_read.sequence) _seq = void;
        bool _rev = void;
    }
}

///
///  Collect per-base information from available tags. 
///  Use $(D arg!TagName) to pass a parameter related to a particular tag.
///
///  Example:
///
/// basesWith!"FZ"(arg!"flowOrder"(flow_order), arg!"keySequence"(key_sequence));
///
template basesWith(Tags...) {
    auto basesWith(R, Args...)(R read, Args args) {
        return PerBaseInfo!(R, Tags)(read, args);
    }
}

/// Provides additional property $(D flow_call).
template FZbaseInfo(R) {

    mixin template resultProperties() {
        /// Current flow call
        ReadFlowCall flow_call() @property const {
            return _flow_call;
        }

        private {
            ReadFlowCall _flow_call;
        }
    }

    mixin template rangeMethods() {

        private {
            ReadFlowCallRange!(BamRead.SequenceResult) _flow_calls = void;
            ReadFlowCall _current_flow_call = void;
            ushort _at = void;

            debug {
                string _read_name;
            }
        }

        void setup(Args...)(const ref R read, Args args) 
        {
            string flow_order = void;
            string key_sequence = void;

            debug {
                _read_name = read.read_name.idup;
            }

            enum flowOrderExists = staticIndexOf!(MixinArg!(string, "flowOrder"), Args);
            enum keySequenceExists = staticIndexOf!(MixinArg!(string, "keySequence"), Args);
            static assert(flowOrderExists != -1, `Flow order must be provided via arg!"flowOrder"`);
            static assert(keySequenceExists != -1, `Flow order must be provided via arg!"keySequence"`);

            foreach (arg; args) {
                static if(is(typeof(arg) == MixinArg!(string, "flowOrder")))
                    flow_order = arg;

                static if(is(typeof(arg) == MixinArg!(string, "keySequence")))
                    key_sequence = arg;
            }

            _at = 0;

            _flow_calls = readFlowCalls(read, flow_order, key_sequence);
            if (!_flow_calls.empty) {
                _current_flow_call = _flow_calls.front;
            }
        }

        void populate(Result)(ref Result result) {
            result._flow_call = _current_flow_call;

            debug {
                if ((_rev && result.base != result._flow_call.base.complement)
                    || (!_rev && result.base != result._flow_call.base)) {
                    import std.stdio;
                    stderr.writeln("invalid flow call at ", _read_name, ": ", result.position);
                }
            }
        }

        void update(const ref R read) 
        {
            ++_at;
            if (_at == _current_flow_call.length) {
                _flow_calls.popFront();
                if (!_flow_calls.empty) {
                    _current_flow_call = _flow_calls.front;
                    _at = 0;
                }
            }
        }

        void copy(Range)(ref Range source, ref Range target) {
            target.FZ._flow_calls = source._flow_calls.save();
            target.FZ._at = source.FZ._at;
            target.FZ._current_flow_call = source._current_flow_call;

            debug {
                target._read_name = _read_name;
            }
        }
    }
}

/// Provides additional properties
///     * position
///     * cigar_operation
template CIGARbaseInfo(R) {

    mixin template resultProperties() {
        /// Current CIGAR operation
        CigarOperation cigar_operation() @property const {
            return _cigar_operation;
        }

        /// Position of the corresponding base on the reference.
        /// If current CIGAR operation is not one of 'M', '=', 'X',
        /// returns the position of the previous valid base.
        ulong position() @property const {
            return _reference_position;
        }

        private {
            CigarOperation _cigar_operation;
            ulong _reference_position;
        }
    }

    mixin template rangeMethods() {

        private {
            const(CigarOperation)[] _cigar = void;
            long _index = void;
            CigarOperation _current_cigar_op = void;
            ulong _at = void;
            ulong _ref_pos = void;
        }

        void setup(Args...)(const ref R read, Args) 
        {
            _cigar = read.cigar;

            _index = _rev ? _cigar.length : -1;
            _ref_pos = _rev ? (read.position + read.basesCovered() - 1)
                            : read.position;

            _moveToNextCigarOperator();
        }

        void populate(Result)(ref Result result) {
            result._cigar_operation = _current_cigar_op;
            result._reference_position = _ref_pos;
        }

        void update(const ref R read) 
        {
           ++_at;
           if (_current_cigar_op.is_reference_consuming) {
               _ref_pos += _rev ? -1 : 1;
           }

           if (_at == _current_cigar_op.length) {
               _moveToNextCigarOperator();
           }
        }

        void copy(Range)(const ref Range source, ref Range target) {
            target.CIGAR._cigar = source.CIGAR._cigar;
            target.CIGAR._index = source.CIGAR._index;
            target.CIGAR._current_cigar_op = source.CIGAR._current_cigar_op;
            target.CIGAR._at = source.CIGAR._at;
            target.CIGAR._ref_pos = source.CIGAR._ref_pos;
        }

        private void _moveToNextCigarOperator() {
            _at = 0;
            if (!_rev) {
                for (++_index; _index < _cigar.length; ++_index)
                {
                    _current_cigar_op = _cigar[_index];
                    if (_current_cigar_op.is_query_consuming)
                        break;
                    if (_current_cigar_op.is_reference_consuming)
                        _ref_pos += _current_cigar_op.length;
                }
            } else {
                for (--_index; _index >= 0; --_index) {
                    _current_cigar_op = _cigar[_index];
                    if (_current_cigar_op.is_query_consuming)
                        break;
                    if (_current_cigar_op.is_reference_consuming)
                        _ref_pos += _current_cigar_op.length;
                }
            }
        }
    }
}
