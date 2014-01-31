/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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

module bio.bam.region;

///
struct BamRegion {
    uint ref_id; /// Reference ID in the BAM file
    uint start;  /// 0-based leftmost coordinate (included)
    uint end;    /// 0-based rightmost coordinate (excluded)

    int opCmp(const ref BamRegion other) const nothrow {
        if (this.ref_id > other.ref_id) { return  1; }
        if (this.ref_id < other.ref_id) { return  -1; }
	
        if (this.start > other.start) { return  1; }
        if (this.start < other.start) { return  -1; }

	if (this.end > other.end) { return  1; }
        if (this.end < other.end) { return  -1; }

        return 0;
    }
}
