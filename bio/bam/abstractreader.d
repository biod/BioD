/*
    This file is part of BioD.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>

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
module bio.bam.abstractreader;

import bio.sam.header;
import bio.bam.read;
import bio.bam.referenceinfo;

public import std.range;

/// Common interface for $(DPREF2 bam, reader, BamReader) 
/// and $(DPREF2 sam, reader, SamReader).
interface IBamSamReader {
    /// SAM header
    bio.sam.header.SamHeader header() @property;

    /// Information about reference sequences
    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const;

    /// All reads in the file
    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property;

    /// Filename
    string filename() @property const;
} 
