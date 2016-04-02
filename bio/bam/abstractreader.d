/*
    This file is part of BioD.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

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
    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const nothrow;

    /// All reads in the file
    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property;

    /// Filename
    string filename() @property const;

    ///
    void assumeSequentialProcessing();
} 
