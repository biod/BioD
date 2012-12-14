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
import bio.core.sequence;

import bio.bam.read;
import bio.bam.tagvalue;
import bio.bam.fz.flowcall;
import bio.bam.md.core;

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
    //          Current base of the result is updated before the call.
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
        /// Actual read base, with strand taken into account.
        Base base;
        alias base this;

        string opCast(T)() if (is(T == string))
        {
            return to!string(base);
        }

        bool opEquals(T)(T base) if (is(Unqual!T == Base))
        {
            return this.base == base;
        }

        bool opEquals(T)(T result) if (is(Unqual!T == Result))
        {
            return this == result;
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
        _rev = read.is_reverse_strand;
        _seq = reversableRange!complementBase(read.sequence, _rev);

        foreach (t; Extensions) {
            setup!t(read, args);
        }
    }

    bool empty() @property {
        return _seq.empty;
    }

    Result front() @property {
        Result r = void;
        r.base = _seq.front;
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
        r._seq = _seq.save;
        r._rev = _rev;
        foreach (t; Extensions)
            copy!t(r);
        return r;
    }

    private void moveToNextBase() {

        foreach (t; Extensions) {
            update!t();
        }

        _seq.popFront();
    }

    /// Returns true if the read is reverse strand,
    /// and false otherwise.
    bool reverse_strand() @property const {
        return _rev;
    }

    private {
        bool _rev = void;
        R _read = void;
        ReversableRange!(complementBase, typeof(_read.sequence)) _seq = void;
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

/// Provides additional properties
///     * reference_base
///     * md_operation
template MDbaseInfo(R) {

    mixin template resultProperties() {
        /// If current CIGAR operation is reference consuming,
        /// returns reference base at this position, otherwise
        /// returns '-'.
        ///
        /// If read is on '-' strand, the result will be
        /// complementary base.
        char reference_base() @property {
            return _ref_base;
        }

        /// Current MD operation
        MdOperation md_operation() @property {
            return _md_op;
        }

        private {
            char _ref_base = void;
            MdOperation _md_op = void;
        }
    }

    mixin template rangeMethods() {

        private {
            ReversableRange!(reverseMdOp, MdOperationRange) _md_ops = void;
            uint _match; // remaining length of current match operation
            MdOperation _md_front = void;
        }

        void setup(Args...)(const ref R read, Args args)
        {
            auto md = read["MD"];
            auto md_str = *(cast(string*)&md);
            _md_ops = reversableRange!reverseMdOp(mdOperations(md_str),
                                                  read.is_reverse_strand);
          
            while (!_md_ops.empty)
            {
                _md_front = _md_ops.front;
                _md_ops.popFront();
                if (!_md_front.is_deletion) {
                    if (_md_front.is_match) {
                        _match = _md_front.match;
                    }
                    break;
                }
            }
        }

        void populate(Result)(ref Result result)
        {
            if (!current_cigar_operation.is_reference_consuming)
            {
                result._ref_base = '-';
                return;
            }

            MdOperation op = _md_front;
            if (op.is_mismatch)
                result._ref_base = op.mismatch.asCharacter;
            else if (op.is_match) {
                result._ref_base = result.base.asCharacter;
            }
            else assert(0);

            result._md_op = op;
        }

        void update(const ref R read)
        {
            if (!current_cigar_operation.is_reference_consuming)
                return;

            if (_md_front.is_mismatch)
            {
                if (_md_ops.empty)
                    return;

                _md_front = _md_ops.front;    
                _md_ops.popFront();
            }
            else if (_md_front.is_match)
            {
                --_match;
                if (_match == 0 && !_md_ops.empty) {
                    _md_front = _md_ops.front;
                    _md_ops.popFront();
                }
            }
            else assert(0);

            while (_md_front.is_deletion) {
                if (_md_ops.empty)
                    return;

                _md_front = _md_ops.front;
                _md_ops.popFront();
            }

            if (_match == 0 && _md_front.is_match)
                _match = _md_front.match;
        }

        void copy(Range)(ref Range source, ref Range target)
        {
            target.MD._md_ops = source.MD._md_ops.save;
            target.MD._md_front = source.MD._md_front;
        }
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
                _read_name = read.name.idup;
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
                if (result.base != result._flow_call.base) {
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
///     * cigar_operation_offset
template CIGARbaseInfo(R) {

    mixin template resultProperties() {

        version(CigarExtra)
        {
            /// Current CIGAR operation
            CigarOperation cigar_operation() @property {
                return _cigar[_operation_index];
            }

            /// CIGAR operations before current one
            auto cigar_before() @property {
                return _cigar[0 .. _operation_index];
            }

            /// CIGAR operations after current one
            auto cigar_after() @property {
                return _cigar[_operation_index + 1 .. _cigar.length];
            }
        }
        else
        {
            /// Current CIGAR operation
            CigarOperation cigar_operation() @property const {
                return _current_cigar_op;
            }
        }

        /// Position of the corresponding base on the reference.
        /// If current CIGAR operation is not one of 'M', '=', 'X',
        /// returns the position of the previous mapped base.
        uint position() @property const {
            return _reference_position;
        }

        /// Offset in current CIGAR operation, starting from 0.
        uint cigar_operation_offset() @property const {
            return _cigar_operation_offset;
        }

        private {
            int _operation_index = void;
            uint _reference_position = void;
            uint _cigar_operation_offset = void;
            version (CigarExtra)
            {
                ReversableRange!(identity, const(CigarOperation)[]) _cigar = void;
            }
            else
            {
                CigarOperation _current_cigar_op;
            }
        }
    }

    mixin template rangeMethods() {

        private {
            CigarOperation _current_cigar_op = void;
            int _index = void;
            uint _at = void;
            uint _ref_pos = void;
            ReversableRange!(identity, const(CigarOperation)[]) _cigar = void;
        }

        /// Current CIGAR operation, available to all extensions
        const(CigarOperation) current_cigar_operation() @property const {
            return _current_cigar_op;
        }

        void setup(Args...)(const ref R read, Args) 
        {
            _cigar = reversableRange(read.cigar, read.is_reverse_strand);

            _index = -1;
            _ref_pos = reverse_strand ? (read.position + read.basesCovered() - 1)
                                      : read.position;

            _moveToNextCigarOperator();
            assert(_index >= 0);
        }

        void populate(Result)(ref Result result) {
            result._reference_position = _ref_pos;
            result._cigar_operation_offset = _at;
            version (CigarExtra)
            {
                result._cigar = _cigar;
                result._operation_index = _index;
            }
            else
            {
                result._current_cigar_op = _current_cigar_op;
            }
        }

        void update(const ref R read) 
        {
           ++_at;
           if (_current_cigar_op.is_reference_consuming) {
               _ref_pos += reverse_strand ? -1 : 1;
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
            for (++_index; _index < _cigar.length; ++_index)
            {
                _current_cigar_op = _cigar[_index];
                if (_current_cigar_op.is_query_consuming)
                    break;
                if (_current_cigar_op.is_reference_consuming)
                {
                    if (reverse_strand)
                        _ref_pos -= _current_cigar_op.length;
                    else
                        _ref_pos += _current_cigar_op.length;
                }
            }
        }
    }
}
