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
module bio.maf.block;

///
struct MafSequence {
    ///
    size_t size;
    ///
    size_t start;

    ///
    string source;
    ///
    char strand;

    ///
    size_t source_size;
    ///
    string text;
   
    ///
    char left_status = '.';
    ///
    size_t left_count;

    ///
    char right_status = '.';
    ///
    size_t right_count;

    ///
    bool empty() @property {
        return empty_status != '.';
    }

    ///
    char empty_status = '.';

    ///
    string quality;
}

///
struct MafBlock {
    ///
    MafSequence[] sequences;
    ///
    double score = double.nan;
    ///
    uint pass = 0;
}
