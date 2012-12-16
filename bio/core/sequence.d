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
module bio.core.sequence;

import bio.core.base;

import std.algorithm;
import std.range;
import std.conv;
import std.traits;
public import std.array;

/// Identity function
T identity(T)(auto ref T t) { return t; }

/// Range that allows to unify operations in forward and reverse directions
/// without virtual function call overhead introduced by $(D inputRangeObject).
///
/// $(D reverseTransform) is a function that will be applied to elements
/// if range is iterated backwards.
struct ReversableRange(alias reverseTransform=identity, R)
    if(isBidirectionalRange!R)
{
    private 
    { 
        bool _rev = void;
        R _range = void;
    }

    /// Construct reversable range.
    ///
    /// Params:
    ///     range   =  bidirectional range
    ///     reverse =  if true, all operations on the range will be as if
    ///                $(D retro(range)) was used instead of $(D range).
    this(R range, bool reverse=false) 
    {
        _rev = reverse;
        _range = range;
    }

    /// Bidirectional range primitives
    bool empty() @property
    {
        return _range.empty;
    }

    /// ditto
    auto front() @property
    {
        return _rev ? reverseTransform(_range.back) : _range.front;
    }

    /// ditto
    auto back() @property
    {
        return _rev ? reverseTransform(_range.front) : _range.back;
    }

    /// ditto
    void popFront() 
    { 
        if (_rev) 
            _range.popBack(); 
        else 
            _range.popFront(); 
    }

    /// ditto
    void popBack()
    {
        if (_rev)
            _range.popFront();
        else
            _range.popBack();
    }

    /// ditto
    auto save() @property
    {
        return ReversableRange(_range.save, _rev);
    }

    /// Reverse of this range
    ReversableRange reverse() @property {
        return ReversableRange(_range.save, !_rev);
    }

    static if(hasLength!R) 
    {
        /// If source range has length, the result also has length
        size_t length() @property
        {
            return _range.length;
        }
    }

    static if(isRandomAccessRange!R)
    {
        /// If source range is a random access range, $(D opIndex) is defined
        auto opIndex(size_t index)
        {
            if (_rev)
                return reverseTransform(_range[_range.length - 1 - index]);
            else
                return _range[index];
        }
    }

    static if(hasSlicing!R)
    {
        /// Slicing is also propagated
        auto opSlice(size_t from, size_t to)
        {
            if (_rev)
            {
                auto len = _range.length;
                //
                //  [b, e) -> (len - 1 - e, len - 1 - b] ~ [len - e, len - b)
                //
                return ReversableRange(_range[len - to .. len - from], true);
            }
            else
                return ReversableRange(_range[from .. to], false);
        }
    }
}

/// Create reversable range from bidirectional one.
ReversableRange!(reverseTransform, R)
reversableRange(alias reverseTransform=identity, R)(R range, bool reverse=false)
{
    return typeof(return)(range, reverse);
}

unittest {

    auto bidir_range = [1, 2, 3, 4, 5];
    auto rev = reversableRange(bidir_range[], true);

    assert(rev.front == 5);
    assert(rev[2] == 3);
    rev.popFront();
    assert(rev.back == 1);
    assert(rev.front == 4);
    assert(equal(rev[1 .. 3], [3, 2]));

    // Here. That's the whole point.
    // One can't do the same with $(D retro)
    // without using $(D inputRangeObject),
    // but that kills performance because
    // virtual calls can not be inlined.
    rev = reversableRange(bidir_range[], false);

    assert(rev.front == 1);
    assert(equal(rev[1 .. 3], [2, 3]));
}

/// Sequence of bases. Element of reversed range will be complemented.
template Sequence(R)
{
    alias ReversableRange!(complementBase, R) Sequence;
}

/// Returns an object very similar to string, but sliceable.
/// Tricks std.traits.isNarrowString.
auto sliceableString(string s) {
    return map!"cast(char)a"(cast(ubyte[])s);
}

///
alias ReturnType!sliceableString SliceableString;

/// Create nucleotide sequence from bidirectional base range.
auto nucleotideSequence(R)(R bases, bool reverse=false)
    if(isBidirectionalRange!R)
{

    static if(isNarrowString!R)
    {
        return nucleotideSequence(sliceableString(bases), reverse);
    } 
    else static if(is(Unqual!(ElementType!R) == char) || 
              is(Unqual!(ElementType!R) == dchar))
    {
        return nucleotideSequence(map!(charToBase!Base16)(bases), reverse);
    }
    else
    {
        return Sequence!R(bases, reverse);
    }
}

///
alias ReturnType!(nucleotideSequence!SliceableString) NucleotideSequence;

unittest {
    auto seq0 = nucleotideSequence("ACGTACGT");

    // reverse-complement
    assert(equal(seq0.reverse[2 .. 6], "GTAC"));

    auto seq1 = nucleotideSequence(seq0, true);
    assert(equal(seq1[1 .. 5], "CGTA"));
    assert(equal(seq1, map!complementBase(retro(seq0))));

    seq1 = nucleotideSequence(seq0, false);
    assert(equal(seq1, seq0));
}
