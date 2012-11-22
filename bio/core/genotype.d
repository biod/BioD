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
module bio.core.genotype;

import bio.core.base;
import bio.core.tinymap;

/// Holds ordered pair of two alleles
struct DiploidGenotype(B) {

    mixin TinyMapInterface!(B.ValueSetSize ^^ 2);

    private static ubyte _getCode(B b1, B b2) {
        auto c1 = b1.internal_code;
        auto c2 = b2.internal_code;
        return cast(ubyte)(c1 * B.ValueSetSize + c2);
    }

    /// Construct a genotype from two bases
    /// Every ambiguous base gets converted to 'N' internally.
    this(B b1, B b2) {
        _code = _getCode(b1, b2);
    }

    /// Construct homozygous genotype
    this(B b) {
        _code = _getCode(b, b);
    }

    /// First allele
    B base1() @property const {
        return B.fromInternalCode(_code / B.ValueSetSize);
    }

    /// Second allele
    B base2() @property const {
        return B.fromInternalCode(_code % B.ValueSetSize);
    }

    ///
    bool is_heterozygous() @property const {
        return base1 != base2;
    }

    ///
    bool is_homozygous() @property const {
        return base1 == base2;
    }

    /// String representation B1|B2 (TODO: add phasing in future?)
    string toString() const {
        return base1 ~ "|" ~ base2;
    }
}

/// Create an instance of DiploidGenotype
auto diploidGenotype(B...)(B bases) {
    return DiploidGenotype!(B[0])(bases);
}

unittest {
    auto g1 = diploidGenotype(Base('C'), Base('W'));
    assert(g1.base1 == 'C');
    assert(g1.base2 == 'W');

    // By default, Base5 is used
    auto g2 = diploidGenotype(Base5('C'));
    assert(g2.base1 == g2.base2);

    // Both bases must be of the same type
    static assert(!__traits(compiles, diploidGenotype(Base5('T'), Base16('D'))));
}
