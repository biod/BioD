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
module bio.bam.referenceinfo;

import std.stream;
import std.exception;
import std.array;

/**
  Stores basic information about reference sequence.
 */
struct ReferenceSequenceInfo {
    private {
        string _name;
        int _length;
    }

    /// Reference sequence name
	/// (null byte is guaranteed to follow the returned slice)
    string name() @property const {
        return _name[0 .. $ - 1];
    }
   
    /// Reference sequence length
    int length() @property const {
        return _length;
    }

    ///
    this(string name, int length) {
        _name = name ~ '\0';
        _length = length;
    }

    /// Constructs the structure from input stream
    this(ref Stream stream) {
        int l_name; // length of the reference name plus one
        stream.read(l_name);
        _name = stream.readString(l_name).idup; // keep '\0' at the end
        stream.read(_length);
    }
}
