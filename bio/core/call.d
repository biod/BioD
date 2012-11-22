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
module bio.core.call;

import bio.core.base;
import bio.core.genotype;

/// A genotype call
struct Call(alias Gt, B) 
{
    alias Gt!B G;

    private {
        string _sample = void;
        string _chr = void;
        ulong _pos = void;
        B _refbase = void;
        G _gt = void;
        float _qual = void;
    }

    /// Constructor
    this(string sample, string chr, ulong pos,
         B refbase, G genotype, float quality=float.nan)
    {
        _sample = sample;
        _chr = chr;
        _pos = pos;
        _refbase = refbase;
        _gt = genotype;
        _qual = quality;
    }

    /// Sample name
    string sample() @property const {
        return _sample;
    }

    /// Chromosome name
    string chromosome() @property const {
        return _chr;
    }

    /// 0-based position on the reference
    ulong position() @property const {
        return _pos;
    }

    /// Reference base at the site
    B reference_base() @property const {
        return _refbase;
    }

    /// Most probable genotype
    ref const(G) genotype() @property const {
        return _gt;
    }

    ///
    alias genotype this; 

    /// Phred-scaled quality score. If unknown, set to NaN.
    float quality() @property const {
        return _qual;
    }

    /// Returns true if this call is not a reference one.
    bool is_variant() @property const {
        return _gt != G(_refbase);
    }
}

alias Call!(DiploidGenotype, Base5) DiploidCall5;
alias Call!(DiploidGenotype, Base16) DiploidCall;
alias DiploidCall DiploidCall16;

unittest {
    auto call = DiploidCall("NA01234", "chr10", 543210,
                            Base('T'), diploidGenotype(Base('C'), Base('T')),
                            47.0);

    assert(call.is_variant);
    assert(call.is_heterozygous);
    assert(call.reference_base == 'T');
}
