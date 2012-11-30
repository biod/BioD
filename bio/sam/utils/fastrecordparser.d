module bio.sam.utils.fastrecordparser;

#line 1 "sam_alignment.rl"
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

#line 24 "sam_alignment.d"
static const int sam_alignment_start = 0;
static const int sam_alignment_first_final = 135;
static const int sam_alignment_error = -1;

static const int sam_alignment_en_alignment = 0;


#line 234 "sam_alignment.rl"


import bio.sam.header;
import bio.bam.read;
import bio.bam.tagvalue;
import bio.bam.utils.tagstoragebuilder;

import std.array;
import std.conv;
import std.typecons;
import std.outbuffer;
import std.c.stdlib;

class AlignmentBuildStorage {
    Appender!(CigarOperation[]) cigar_appender;
    OutBuffer outbuffer;
    TagStorageBuilder tag_storage_builder;

    this() {
        cigar_appender = appender!(CigarOperation[])();
        outbuffer = new OutBuffer();
        tag_storage_builder = TagStorageBuilder.create();
    }

    void clear() {
        cigar_appender.clear();
        tag_storage_builder.clear();
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
    }
}

BamRead parseAlignmentLine(string line, SamHeader header, 
                             AlignmentBuildStorage b=null) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    if (b is null) {
        b = new AlignmentBuildStorage();
    } else {
        b.clear();
    }

    byte current_sign = 1;

    size_t read_name_beg; // position of beginning of QNAME
    size_t read_name_end; // position past the end of QNAME

    size_t sequence_beg; // position of SEQ start
    string sequence;     // SEQ

    uint cigar_op_len;   // length of CIGAR operation
    char cigar_op_chr;   // CIGAR operation

    size_t cigar_op_len_start; // position of start of CIGAR operation
    
    auto cigar = b.cigar_appender;

    long int_value;                      // for storing temporary integers
    float float_value;                   // for storing temporary floats
    size_t float_beg;                    // position of start of current float
    auto outbuffer = b.outbuffer;        // used to build tag values which hold arrays
    char arraytype;                      // type of last array tag value

    ushort flag;
    uint pos;
    uint mate_pos;
    ubyte mapping_quality; 
    int template_length;
    ubyte* qual_ptr = null;
    size_t qual_index; 

    string current_tag;
    Value current_tagvalue;

    size_t tag_key_beg, tagvalue_beg;
    size_t rname_beg, rnext_beg;

    int ref_id = -1;
    int mate_ref_id = -1;
    
    auto builder = b.tag_storage_builder;

    
#line 119 "sam_alignment.d"
	{
	cs = sam_alignment_start;
	}

#line 320 "sam_alignment.rl"
    
#line 126 "sam_alignment.d"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
goto case; case 0:
	if ( (*p) == 9u )
		goto tr1;
	if ( (*p) > 63u ) {
		if ( 65u <= (*p) && (*p) <= 126u )
			goto tr2;
	} else if ( (*p) >= 33u )
		goto tr2;
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
goto case; case 1:
	if ( (*p) == 9u )
		goto st2;
	goto st1;
tr1:
#line 43 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
#line 44 "sam_alignment.rl"
	{ read_name_end = p - line.ptr; }
	goto st2;
tr150:
#line 44 "sam_alignment.rl"
	{ read_name_end = p - line.ptr; }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
goto case; case 2:
#line 162 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st4;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr6;
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
goto case; case 3:
	if ( (*p) == 9u )
		goto st4;
	goto st3;
tr132:
#line 48 "sam_alignment.rl"
	{ flag = to!ushort(int_value); }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
goto case; case 4:
#line 183 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st6;
	if ( (*p) < 43u ) {
		if ( 33u <= (*p) && (*p) <= 41u )
			goto tr9;
	} else if ( (*p) > 60u ) {
		if ( 62u <= (*p) && (*p) <= 126u )
			goto tr9;
	} else
		goto tr9;
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
goto case; case 5:
	if ( (*p) == 9u )
		goto st6;
	goto st5;
tr130:
#line 55 "sam_alignment.rl"
	{
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]); 
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
goto case; case 6:
#line 212 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st8;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr12;
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
goto case; case 7:
	if ( (*p) == 9u )
		goto st8;
	goto st7;
tr112:
#line 49 "sam_alignment.rl"
	{ pos = to!uint(int_value); }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
goto case; case 8:
#line 233 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st10;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr15;
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
goto case; case 9:
	if ( (*p) == 9u )
		goto st10;
	goto st9;
tr94:
#line 50 "sam_alignment.rl"
	{ mapping_quality = to!ubyte(int_value); }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
goto case; case 10:
#line 254 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st12;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr18;
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
goto case; case 11:
	if ( (*p) == 9u )
		goto st12;
	goto st11;
tr92:
#line 66 "sam_alignment.rl"
	{ 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
goto case; case 12:
#line 277 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st14;
		case 61u: goto st59;
		default: break;
	}
	if ( (*p) > 41u ) {
		if ( 43u <= (*p) && (*p) <= 126u )
			goto tr21;
	} else if ( (*p) >= 33u )
		goto tr21;
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
goto case; case 13:
	if ( (*p) == 9u )
		goto st14;
	goto st13;
tr71:
#line 78 "sam_alignment.rl"
	{
        mate_ref_id = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
    }
	goto st14;
tr73:
#line 73 "sam_alignment.rl"
	{
        mate_ref_id = ref_id;
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
goto case; case 14:
#line 312 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st16;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr25;
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
goto case; case 15:
	if ( (*p) == 9u )
		goto st16;
	goto st15;
tr53:
#line 85 "sam_alignment.rl"
	{ mate_pos = to!uint(int_value); }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
goto case; case 16:
#line 333 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st18;
		case 43u: goto tr28;
		case 45u: goto tr28;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr29;
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
goto case; case 17:
	if ( (*p) == 9u )
		goto st18;
	goto st17;
tr35:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 86 "sam_alignment.rl"
	{ template_length = to!int(int_value); }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
goto case; case 18:
#line 360 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st135;
		case 46u: goto tr32;
		case 61u: goto tr32;
		default: break;
	}
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto tr32;
	} else if ( (*p) >= 65u )
		goto tr32;
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
goto case; case 19:
	if ( (*p) == 9u )
		goto st135;
	goto st19;
tr33:
#line 92 "sam_alignment.rl"
	{ sequence = line[sequence_beg .. p - line.ptr]; }
	goto st135;
st135:
	if ( ++p == pe )
		goto _test_eof135;
goto case; case 135:
#line 388 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr154;
	goto st136;
st136:
	if ( ++p == pe )
		goto _test_eof136;
goto case; case 136:
	if ( (*p) == 9u )
		goto st137;
	goto st136;
tr166:
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr173:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 180 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr202:
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 193 "sam_alignment.rl"
	{ 
        outbuffer.write(float_value); 
    }
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr212:
#line 161 "sam_alignment.rl"
	{
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]);
        current_tagvalue.setHexadecimalFlag();
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr216:
#line 157 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]); 
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr228:
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 153 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(float_value);
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
tr238:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 127 "sam_alignment.rl"
	{ 
        if (int_value < 0) {
            if (int_value >= byte.min) {
                current_tagvalue = Value(to!byte(int_value));
            } else if (int_value >= short.min) {
                current_tagvalue = Value(to!short(int_value));
            } else if (int_value >= int.min) {
                current_tagvalue = Value(to!int(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                current_tagvalue = Value(to!ubyte(int_value));
            } else if (int_value <= ushort.max) {
                current_tagvalue = Value(to!ushort(int_value));
            } else if (int_value <= uint.max) {
                current_tagvalue = Value(to!uint(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	goto st137;
st137:
	if ( ++p == pe )
		goto _test_eof137;
goto case; case 137:
#line 525 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto tr155;
	} else if ( (*p) >= 65u )
		goto tr155;
	goto st136;
tr155:
#line 223 "sam_alignment.rl"
	{ tag_key_beg = p - line.ptr; }
	goto st138;
st138:
	if ( ++p == pe )
		goto _test_eof138;
goto case; case 138:
#line 542 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto st139;
	} else if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto st139;
	} else
		goto st139;
	goto st136;
st139:
	if ( ++p == pe )
		goto _test_eof139;
goto case; case 139:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto tr157;
		default: break;
	}
	goto st136;
tr157:
#line 224 "sam_alignment.rl"
	{ current_tag = line[tag_key_beg .. p - line.ptr]; }
	goto st140;
st140:
	if ( ++p == pe )
		goto _test_eof140;
goto case; case 140:
#line 572 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 65u: goto st141;
		case 66u: goto st144;
		case 72u: goto st181;
		case 90u: goto st184;
		case 102u: goto st187;
		case 105u: goto st201;
		default: break;
	}
	goto st136;
st141:
	if ( ++p == pe )
		goto _test_eof141;
goto case; case 141:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st142;
		default: break;
	}
	goto st136;
st142:
	if ( ++p == pe )
		goto _test_eof142;
goto case; case 142:
	if ( (*p) == 9u )
		goto st137;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr165;
	goto st136;
tr165:
#line 126 "sam_alignment.rl"
	{ current_tagvalue = Value((*p)); }
	goto st143;
st143:
	if ( ++p == pe )
		goto _test_eof143;
goto case; case 143:
#line 611 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr166;
	goto st136;
st144:
	if ( ++p == pe )
		goto _test_eof144;
goto case; case 144:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st145;
		default: break;
	}
	goto st136;
st145:
	if ( ++p == pe )
		goto _test_eof145;
goto case; case 145:
	switch( (*p) ) {
		case 9u: goto st137;
		case 67u: goto tr168;
		case 73u: goto tr168;
		case 83u: goto tr168;
		case 99u: goto tr168;
		case 102u: goto tr169;
		case 105u: goto tr168;
		case 115u: goto tr168;
		default: break;
	}
	goto st136;
tr168:
#line 170 "sam_alignment.rl"
	{
        // it might be not the best idea to use outbuffer;
        // the better idea might be two-pass approach
        // when first pass is for counting commas, and
        // the second is for filling allocated array
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
        arraytype = (*p);
    }
	goto st146;
st146:
	if ( ++p == pe )
		goto _test_eof146;
goto case; case 146:
#line 657 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 44u: goto st147;
		default: break;
	}
	goto st136;
tr174:
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 180 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }
	goto st147;
st147:
	if ( ++p == pe )
		goto _test_eof147;
goto case; case 147:
#line 685 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto tr171;
		case 45u: goto tr171;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr172;
	goto st136;
tr171:
#line 23 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st148;
st148:
	if ( ++p == pe )
		goto _test_eof148;
goto case; case 148:
#line 703 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr172;
	goto st136;
tr172:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st149;
st149:
	if ( ++p == pe )
		goto _test_eof149;
goto case; case 149:
#line 719 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr175;
	goto st136;
tr175:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st150;
st150:
	if ( ++p == pe )
		goto _test_eof150;
goto case; case 150:
#line 736 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr176;
	goto st136;
tr176:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st151;
st151:
	if ( ++p == pe )
		goto _test_eof151;
goto case; case 151:
#line 753 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr177;
	goto st136;
tr177:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st152;
st152:
	if ( ++p == pe )
		goto _test_eof152;
goto case; case 152:
#line 770 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr178;
	goto st136;
tr178:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st153;
st153:
	if ( ++p == pe )
		goto _test_eof153;
goto case; case 153:
#line 787 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr179;
	goto st136;
tr179:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st154;
st154:
	if ( ++p == pe )
		goto _test_eof154;
goto case; case 154:
#line 804 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr180;
	goto st136;
tr180:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st155;
st155:
	if ( ++p == pe )
		goto _test_eof155;
goto case; case 155:
#line 821 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr181;
	goto st136;
tr181:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st156;
st156:
	if ( ++p == pe )
		goto _test_eof156;
goto case; case 156:
#line 838 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr182;
	goto st136;
tr182:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st157;
st157:
	if ( ++p == pe )
		goto _test_eof157;
goto case; case 157:
#line 855 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr183;
	goto st136;
tr183:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st158;
st158:
	if ( ++p == pe )
		goto _test_eof158;
goto case; case 158:
#line 872 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr184;
	goto st136;
tr184:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st159;
st159:
	if ( ++p == pe )
		goto _test_eof159;
goto case; case 159:
#line 889 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr185;
	goto st136;
tr185:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st160;
st160:
	if ( ++p == pe )
		goto _test_eof160;
goto case; case 160:
#line 906 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr186;
	goto st136;
tr186:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st161;
st161:
	if ( ++p == pe )
		goto _test_eof161;
goto case; case 161:
#line 923 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr187;
	goto st136;
tr187:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st162;
st162:
	if ( ++p == pe )
		goto _test_eof162;
goto case; case 162:
#line 940 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr188;
	goto st136;
tr188:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st163;
st163:
	if ( ++p == pe )
		goto _test_eof163;
goto case; case 163:
#line 957 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr189;
	goto st136;
tr189:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st164;
st164:
	if ( ++p == pe )
		goto _test_eof164;
goto case; case 164:
#line 974 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr190;
	goto st136;
tr190:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st165;
st165:
	if ( ++p == pe )
		goto _test_eof165;
goto case; case 165:
#line 991 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr191;
	goto st136;
tr191:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st166;
st166:
	if ( ++p == pe )
		goto _test_eof166;
goto case; case 166:
#line 1008 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr173;
		case 44u: goto tr174;
		default: break;
	}
	goto st136;
tr169:
#line 170 "sam_alignment.rl"
	{
        // it might be not the best idea to use outbuffer;
        // the better idea might be two-pass approach
        // when first pass is for counting commas, and
        // the second is for filling allocated array
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
        arraytype = (*p);
    }
	goto st167;
st167:
	if ( ++p == pe )
		goto _test_eof167;
goto case; case 167:
#line 1031 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 44u: goto st168;
		default: break;
	}
	goto st136;
tr203:
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 193 "sam_alignment.rl"
	{ 
        outbuffer.write(float_value); 
    }
	goto st168;
st168:
	if ( ++p == pe )
		goto _test_eof168;
goto case; case 168:
#line 1052 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto tr193;
		case 45u: goto tr193;
		case 46u: goto tr194;
		case 105u: goto tr196;
		case 110u: goto tr197;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr195;
	goto st136;
tr193:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st169;
st169:
	if ( ++p == pe )
		goto _test_eof169;
goto case; case 169:
#line 1073 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 46u: goto st170;
		case 105u: goto st176;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st175;
	goto st136;
tr194:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st170;
st170:
	if ( ++p == pe )
		goto _test_eof170;
goto case; case 170:
#line 1091 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st171;
	goto st136;
st171:
	if ( ++p == pe )
		goto _test_eof171;
goto case; case 171:
	switch( (*p) ) {
		case 9u: goto tr202;
		case 44u: goto tr203;
		case 69u: goto st172;
		case 101u: goto st172;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st171;
	goto st136;
st172:
	if ( ++p == pe )
		goto _test_eof172;
goto case; case 172:
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto st173;
		case 45u: goto st173;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st174;
	goto st136;
st173:
	if ( ++p == pe )
		goto _test_eof173;
goto case; case 173:
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st174;
	goto st136;
st174:
	if ( ++p == pe )
		goto _test_eof174;
goto case; case 174:
	switch( (*p) ) {
		case 9u: goto tr202;
		case 44u: goto tr203;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st174;
	goto st136;
tr195:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st175;
st175:
	if ( ++p == pe )
		goto _test_eof175;
goto case; case 175:
#line 1153 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr202;
		case 44u: goto tr203;
		case 46u: goto st170;
		case 69u: goto st172;
		case 101u: goto st172;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st175;
	goto st136;
tr196:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st176;
st176:
	if ( ++p == pe )
		goto _test_eof176;
goto case; case 176:
#line 1173 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 110u: goto st177;
		default: break;
	}
	goto st136;
st177:
	if ( ++p == pe )
		goto _test_eof177;
goto case; case 177:
	switch( (*p) ) {
		case 9u: goto st137;
		case 102u: goto st178;
		default: break;
	}
	goto st136;
st178:
	if ( ++p == pe )
		goto _test_eof178;
goto case; case 178:
	switch( (*p) ) {
		case 9u: goto tr202;
		case 44u: goto tr203;
		default: break;
	}
	goto st136;
tr197:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st179;
st179:
	if ( ++p == pe )
		goto _test_eof179;
goto case; case 179:
#line 1208 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 97u: goto st180;
		default: break;
	}
	goto st136;
st180:
	if ( ++p == pe )
		goto _test_eof180;
goto case; case 180:
	switch( (*p) ) {
		case 9u: goto st137;
		case 110u: goto st178;
		default: break;
	}
	goto st136;
st181:
	if ( ++p == pe )
		goto _test_eof181;
goto case; case 181:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st182;
		default: break;
	}
	goto st136;
st182:
	if ( ++p == pe )
		goto _test_eof182;
goto case; case 182:
	if ( (*p) == 9u )
		goto st137;
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr211;
	} else if ( (*p) > 70u ) {
		if ( 97u <= (*p) && (*p) <= 102u )
			goto tr211;
	} else
		goto tr211;
	goto st136;
tr211:
#line 151 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	goto st183;
st183:
	if ( ++p == pe )
		goto _test_eof183;
goto case; case 183:
#line 1258 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr212;
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto st183;
	} else if ( (*p) > 70u ) {
		if ( 97u <= (*p) && (*p) <= 102u )
			goto st183;
	} else
		goto st183;
	goto st136;
st184:
	if ( ++p == pe )
		goto _test_eof184;
goto case; case 184:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st185;
		default: break;
	}
	goto st136;
st185:
	if ( ++p == pe )
		goto _test_eof185;
goto case; case 185:
	if ( (*p) == 9u )
		goto st137;
	if ( 32u <= (*p) && (*p) <= 126u )
		goto tr215;
	goto st136;
tr215:
#line 151 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	goto st186;
st186:
	if ( ++p == pe )
		goto _test_eof186;
goto case; case 186:
#line 1297 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr216;
	if ( 32u <= (*p) && (*p) <= 126u )
		goto st186;
	goto st136;
st187:
	if ( ++p == pe )
		goto _test_eof187;
goto case; case 187:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st188;
		default: break;
	}
	goto st136;
st188:
	if ( ++p == pe )
		goto _test_eof188;
goto case; case 188:
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto tr219;
		case 45u: goto tr219;
		case 46u: goto tr220;
		case 105u: goto tr222;
		case 110u: goto tr223;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr221;
	goto st136;
tr219:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st189;
st189:
	if ( ++p == pe )
		goto _test_eof189;
goto case; case 189:
#line 1337 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 46u: goto st190;
		case 105u: goto st196;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st195;
	goto st136;
tr220:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st190;
st190:
	if ( ++p == pe )
		goto _test_eof190;
goto case; case 190:
#line 1355 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st191;
	goto st136;
st191:
	if ( ++p == pe )
		goto _test_eof191;
goto case; case 191:
	switch( (*p) ) {
		case 9u: goto tr228;
		case 69u: goto st192;
		case 101u: goto st192;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st191;
	goto st136;
st192:
	if ( ++p == pe )
		goto _test_eof192;
goto case; case 192:
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto st193;
		case 45u: goto st193;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st194;
	goto st136;
st193:
	if ( ++p == pe )
		goto _test_eof193;
goto case; case 193:
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st194;
	goto st136;
st194:
	if ( ++p == pe )
		goto _test_eof194;
goto case; case 194:
	if ( (*p) == 9u )
		goto tr228;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st194;
	goto st136;
tr221:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st195;
st195:
	if ( ++p == pe )
		goto _test_eof195;
goto case; case 195:
#line 1413 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr228;
		case 46u: goto st190;
		case 69u: goto st192;
		case 101u: goto st192;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st195;
	goto st136;
tr222:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st196;
st196:
	if ( ++p == pe )
		goto _test_eof196;
goto case; case 196:
#line 1432 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 110u: goto st197;
		default: break;
	}
	goto st136;
st197:
	if ( ++p == pe )
		goto _test_eof197;
goto case; case 197:
	switch( (*p) ) {
		case 9u: goto st137;
		case 102u: goto st198;
		default: break;
	}
	goto st136;
st198:
	if ( ++p == pe )
		goto _test_eof198;
goto case; case 198:
	if ( (*p) == 9u )
		goto tr228;
	goto st136;
tr223:
#line 33 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st199;
st199:
	if ( ++p == pe )
		goto _test_eof199;
goto case; case 199:
#line 1464 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st137;
		case 97u: goto st200;
		default: break;
	}
	goto st136;
st200:
	if ( ++p == pe )
		goto _test_eof200;
goto case; case 200:
	switch( (*p) ) {
		case 9u: goto st137;
		case 110u: goto st198;
		default: break;
	}
	goto st136;
st201:
	if ( ++p == pe )
		goto _test_eof201;
goto case; case 201:
	switch( (*p) ) {
		case 9u: goto st137;
		case 58u: goto st202;
		default: break;
	}
	goto st136;
st202:
	if ( ++p == pe )
		goto _test_eof202;
goto case; case 202:
	switch( (*p) ) {
		case 9u: goto st137;
		case 43u: goto tr236;
		case 45u: goto tr236;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr237;
	goto st136;
tr236:
#line 23 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st203;
st203:
	if ( ++p == pe )
		goto _test_eof203;
goto case; case 203:
#line 1512 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr237;
	goto st136;
tr237:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st204;
st204:
	if ( ++p == pe )
		goto _test_eof204;
goto case; case 204:
#line 1528 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr239;
	goto st136;
tr239:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st205;
st205:
	if ( ++p == pe )
		goto _test_eof205;
goto case; case 205:
#line 1542 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr240;
	goto st136;
tr240:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st206;
st206:
	if ( ++p == pe )
		goto _test_eof206;
goto case; case 206:
#line 1556 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr241;
	goto st136;
tr241:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st207;
st207:
	if ( ++p == pe )
		goto _test_eof207;
goto case; case 207:
#line 1570 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr242;
	goto st136;
tr242:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st208;
st208:
	if ( ++p == pe )
		goto _test_eof208;
goto case; case 208:
#line 1584 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr243;
	goto st136;
tr243:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st209;
st209:
	if ( ++p == pe )
		goto _test_eof209;
goto case; case 209:
#line 1598 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr244;
	goto st136;
tr244:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st210;
st210:
	if ( ++p == pe )
		goto _test_eof210;
goto case; case 210:
#line 1612 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr245;
	goto st136;
tr245:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st211;
st211:
	if ( ++p == pe )
		goto _test_eof211;
goto case; case 211:
#line 1626 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr246;
	goto st136;
tr246:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st212;
st212:
	if ( ++p == pe )
		goto _test_eof212;
goto case; case 212:
#line 1640 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr247;
	goto st136;
tr247:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st213;
st213:
	if ( ++p == pe )
		goto _test_eof213;
goto case; case 213:
#line 1654 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr248;
	goto st136;
tr248:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st214;
st214:
	if ( ++p == pe )
		goto _test_eof214;
goto case; case 214:
#line 1668 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr249;
	goto st136;
tr249:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st215;
st215:
	if ( ++p == pe )
		goto _test_eof215;
goto case; case 215:
#line 1682 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr250;
	goto st136;
tr250:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st216;
st216:
	if ( ++p == pe )
		goto _test_eof216;
goto case; case 216:
#line 1696 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr251;
	goto st136;
tr251:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st217;
st217:
	if ( ++p == pe )
		goto _test_eof217;
goto case; case 217:
#line 1710 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr252;
	goto st136;
tr252:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st218;
st218:
	if ( ++p == pe )
		goto _test_eof218;
goto case; case 218:
#line 1724 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr253;
	goto st136;
tr253:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st219;
st219:
	if ( ++p == pe )
		goto _test_eof219;
goto case; case 219:
#line 1738 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr254;
	goto st136;
tr254:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st220;
st220:
	if ( ++p == pe )
		goto _test_eof220;
goto case; case 220:
#line 1752 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr255;
	goto st136;
tr255:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st221;
st221:
	if ( ++p == pe )
		goto _test_eof221;
goto case; case 221:
#line 1766 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr238;
	goto st136;
tr154:
#line 97 "sam_alignment.rl"
	{
        if (sequence.length > 1024) {
            qual_ptr = (new ubyte[sequence.length]).ptr;
        } else {
            qual_ptr = cast(ubyte*)alloca(sequence.length);
            if (!qual_ptr) {
                qual_ptr = (new ubyte[sequence.length]).ptr;
            }
        }
        qual_index = 0;
    }
#line 109 "sam_alignment.rl"
	{
        qual_ptr[qual_index++] = cast(ubyte)((*p) - 33);
    }
	goto st222;
tr256:
#line 109 "sam_alignment.rl"
	{
        qual_ptr[qual_index++] = cast(ubyte)((*p) - 33);
    }
	goto st222;
st222:
	if ( ++p == pe )
		goto _test_eof222;
goto case; case 222:
#line 1798 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st137;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr256;
	goto st136;
tr32:
#line 91 "sam_alignment.rl"
	{ sequence_beg = p - line.ptr; }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
goto case; case 20:
#line 1812 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr33;
		case 46u: goto st20;
		case 61u: goto st20;
		default: break;
	}
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto st20;
	} else if ( (*p) >= 65u )
		goto st20;
	goto st19;
tr28:
#line 23 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
goto case; case 21:
#line 1833 "sam_alignment.d"
	if ( (*p) == 9u )
		goto st18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr29;
	goto st17;
tr29:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
goto case; case 22:
#line 1849 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr36;
	goto st17;
tr36:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
goto case; case 23:
#line 1863 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr37;
	goto st17;
tr37:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
goto case; case 24:
#line 1877 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr38;
	goto st17;
tr38:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
goto case; case 25:
#line 1891 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr39;
	goto st17;
tr39:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
goto case; case 26:
#line 1905 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr40;
	goto st17;
tr40:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
goto case; case 27:
#line 1919 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr41;
	goto st17;
tr41:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st28;
st28:
	if ( ++p == pe )
		goto _test_eof28;
goto case; case 28:
#line 1933 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr42;
	goto st17;
tr42:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st29;
st29:
	if ( ++p == pe )
		goto _test_eof29;
goto case; case 29:
#line 1947 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr43;
	goto st17;
tr43:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st30;
st30:
	if ( ++p == pe )
		goto _test_eof30;
goto case; case 30:
#line 1961 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr44;
	goto st17;
tr44:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
goto case; case 31:
#line 1975 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr45;
	goto st17;
tr45:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st32;
st32:
	if ( ++p == pe )
		goto _test_eof32;
goto case; case 32:
#line 1989 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr46;
	goto st17;
tr46:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st33;
st33:
	if ( ++p == pe )
		goto _test_eof33;
goto case; case 33:
#line 2003 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr47;
	goto st17;
tr47:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st34;
st34:
	if ( ++p == pe )
		goto _test_eof34;
goto case; case 34:
#line 2017 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr48;
	goto st17;
tr48:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st35;
st35:
	if ( ++p == pe )
		goto _test_eof35;
goto case; case 35:
#line 2031 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr49;
	goto st17;
tr49:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
goto case; case 36:
#line 2045 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr50;
	goto st17;
tr50:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
goto case; case 37:
#line 2059 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr51;
	goto st17;
tr51:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
goto case; case 38:
#line 2073 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr52;
	goto st17;
tr52:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st39;
st39:
	if ( ++p == pe )
		goto _test_eof39;
goto case; case 39:
#line 2087 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr35;
	goto st17;
tr25:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
goto case; case 40:
#line 2101 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr54;
	goto st15;
tr54:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
goto case; case 41:
#line 2115 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr55;
	goto st15;
tr55:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st42;
st42:
	if ( ++p == pe )
		goto _test_eof42;
goto case; case 42:
#line 2129 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr56;
	goto st15;
tr56:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st43;
st43:
	if ( ++p == pe )
		goto _test_eof43;
goto case; case 43:
#line 2143 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr57;
	goto st15;
tr57:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
goto case; case 44:
#line 2157 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr58;
	goto st15;
tr58:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
goto case; case 45:
#line 2171 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr59;
	goto st15;
tr59:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st46;
st46:
	if ( ++p == pe )
		goto _test_eof46;
goto case; case 46:
#line 2185 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr60;
	goto st15;
tr60:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st47;
st47:
	if ( ++p == pe )
		goto _test_eof47;
goto case; case 47:
#line 2199 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr61;
	goto st15;
tr61:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
goto case; case 48:
#line 2213 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr62;
	goto st15;
tr62:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
goto case; case 49:
#line 2227 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr63;
	goto st15;
tr63:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
goto case; case 50:
#line 2241 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr64;
	goto st15;
tr64:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
goto case; case 51:
#line 2255 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr65;
	goto st15;
tr65:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st52;
st52:
	if ( ++p == pe )
		goto _test_eof52;
goto case; case 52:
#line 2269 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr66;
	goto st15;
tr66:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
goto case; case 53:
#line 2283 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr67;
	goto st15;
tr67:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
goto case; case 54:
#line 2297 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr68;
	goto st15;
tr68:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st55;
st55:
	if ( ++p == pe )
		goto _test_eof55;
goto case; case 55:
#line 2311 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr69;
	goto st15;
tr69:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
goto case; case 56:
#line 2325 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr70;
	goto st15;
tr70:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
goto case; case 57:
#line 2339 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr53;
	goto st15;
tr21:
#line 77 "sam_alignment.rl"
	{ rnext_beg = p - line.ptr; }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
goto case; case 58:
#line 2351 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr71;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st58;
	goto st13;
st59:
	if ( ++p == pe )
		goto _test_eof59;
goto case; case 59:
	if ( (*p) == 9u )
		goto tr73;
	goto st13;
tr18:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st60;
tr93:
#line 66 "sam_alignment.rl"
	{ 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
goto case; case 60:
#line 2384 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr74;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr74:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
goto case; case 61:
#line 2411 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr76;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr76:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
goto case; case 62:
#line 2438 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr77;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr77:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
goto case; case 63:
#line 2465 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr78;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr78:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
goto case; case 64:
#line 2492 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr79;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr79:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
goto case; case 65:
#line 2519 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr80;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr80:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
goto case; case 66:
#line 2546 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr81;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr81:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
goto case; case 67:
#line 2573 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr82;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr82:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
goto case; case 68:
#line 2600 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr83;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr83:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
goto case; case 69:
#line 2627 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr84;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr84:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
goto case; case 70:
#line 2654 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr85;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr85:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
goto case; case 71:
#line 2681 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr86;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr86:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st72;
st72:
	if ( ++p == pe )
		goto _test_eof72;
goto case; case 72:
#line 2708 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr87;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr87:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
goto case; case 73:
#line 2735 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr88;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr88:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
goto case; case 74:
#line 2762 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr89;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr89:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st75;
st75:
	if ( ++p == pe )
		goto _test_eof75;
goto case; case 75:
#line 2789 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr90;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr90:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st76;
st76:
	if ( ++p == pe )
		goto _test_eof76;
goto case; case 76:
#line 2816 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr91;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else
		goto tr75;
	goto st11;
tr91:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st77;
st77:
	if ( ++p == pe )
		goto _test_eof77;
goto case; case 77:
#line 2843 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto st12;
		case 61u: goto tr75;
		case 68u: goto tr75;
		case 80u: goto tr75;
		case 83u: goto tr75;
		case 88u: goto tr75;
		default: break;
	}
	if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr75;
	} else if ( (*p) >= 72u )
		goto tr75;
	goto st11;
tr75:
#line 64 "sam_alignment.rl"
	{ cigar_op_len = to!uint(int_value); }
#line 65 "sam_alignment.rl"
	{ cigar_op_chr = (*p); }
	goto st78;
st78:
	if ( ++p == pe )
		goto _test_eof78;
goto case; case 78:
#line 2869 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr92;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr93;
	goto st11;
tr15:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
goto case; case 79:
#line 2885 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr95;
	goto st9;
tr95:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st80;
st80:
	if ( ++p == pe )
		goto _test_eof80;
goto case; case 80:
#line 2899 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr96;
	goto st9;
tr96:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st81;
st81:
	if ( ++p == pe )
		goto _test_eof81;
goto case; case 81:
#line 2913 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr97;
	goto st9;
tr97:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
goto case; case 82:
#line 2927 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr98;
	goto st9;
tr98:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
goto case; case 83:
#line 2941 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr99;
	goto st9;
tr99:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st84;
st84:
	if ( ++p == pe )
		goto _test_eof84;
goto case; case 84:
#line 2955 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr100;
	goto st9;
tr100:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
goto case; case 85:
#line 2969 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr101;
	goto st9;
tr101:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st86;
st86:
	if ( ++p == pe )
		goto _test_eof86;
goto case; case 86:
#line 2983 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr102;
	goto st9;
tr102:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
goto case; case 87:
#line 2997 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr103;
	goto st9;
tr103:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
goto case; case 88:
#line 3011 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr104;
	goto st9;
tr104:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st89;
st89:
	if ( ++p == pe )
		goto _test_eof89;
goto case; case 89:
#line 3025 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr105;
	goto st9;
tr105:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st90;
st90:
	if ( ++p == pe )
		goto _test_eof90;
goto case; case 90:
#line 3039 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr106;
	goto st9;
tr106:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st91;
st91:
	if ( ++p == pe )
		goto _test_eof91;
goto case; case 91:
#line 3053 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr107;
	goto st9;
tr107:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st92;
st92:
	if ( ++p == pe )
		goto _test_eof92;
goto case; case 92:
#line 3067 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr108;
	goto st9;
tr108:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st93;
st93:
	if ( ++p == pe )
		goto _test_eof93;
goto case; case 93:
#line 3081 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr109;
	goto st9;
tr109:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st94;
st94:
	if ( ++p == pe )
		goto _test_eof94;
goto case; case 94:
#line 3095 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr110;
	goto st9;
tr110:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st95;
st95:
	if ( ++p == pe )
		goto _test_eof95;
goto case; case 95:
#line 3109 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr111;
	goto st9;
tr111:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st96;
st96:
	if ( ++p == pe )
		goto _test_eof96;
goto case; case 96:
#line 3123 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr94;
	goto st9;
tr12:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st97;
st97:
	if ( ++p == pe )
		goto _test_eof97;
goto case; case 97:
#line 3137 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr113;
	goto st7;
tr113:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st98;
st98:
	if ( ++p == pe )
		goto _test_eof98;
goto case; case 98:
#line 3151 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr114;
	goto st7;
tr114:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st99;
st99:
	if ( ++p == pe )
		goto _test_eof99;
goto case; case 99:
#line 3165 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr115;
	goto st7;
tr115:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st100;
st100:
	if ( ++p == pe )
		goto _test_eof100;
goto case; case 100:
#line 3179 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr116;
	goto st7;
tr116:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st101;
st101:
	if ( ++p == pe )
		goto _test_eof101;
goto case; case 101:
#line 3193 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr117;
	goto st7;
tr117:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st102;
st102:
	if ( ++p == pe )
		goto _test_eof102;
goto case; case 102:
#line 3207 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr118;
	goto st7;
tr118:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st103;
st103:
	if ( ++p == pe )
		goto _test_eof103;
goto case; case 103:
#line 3221 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr119;
	goto st7;
tr119:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st104;
st104:
	if ( ++p == pe )
		goto _test_eof104;
goto case; case 104:
#line 3235 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr120;
	goto st7;
tr120:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st105;
st105:
	if ( ++p == pe )
		goto _test_eof105;
goto case; case 105:
#line 3249 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr121;
	goto st7;
tr121:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
goto case; case 106:
#line 3263 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr122;
	goto st7;
tr122:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
goto case; case 107:
#line 3277 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr123;
	goto st7;
tr123:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st108;
st108:
	if ( ++p == pe )
		goto _test_eof108;
goto case; case 108:
#line 3291 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr124;
	goto st7;
tr124:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st109;
st109:
	if ( ++p == pe )
		goto _test_eof109;
goto case; case 109:
#line 3305 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr125;
	goto st7;
tr125:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st110;
st110:
	if ( ++p == pe )
		goto _test_eof110;
goto case; case 110:
#line 3319 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr126;
	goto st7;
tr126:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st111;
st111:
	if ( ++p == pe )
		goto _test_eof111;
goto case; case 111:
#line 3333 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr127;
	goto st7;
tr127:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st112;
st112:
	if ( ++p == pe )
		goto _test_eof112;
goto case; case 112:
#line 3347 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr128;
	goto st7;
tr128:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st113;
st113:
	if ( ++p == pe )
		goto _test_eof113;
goto case; case 113:
#line 3361 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr129;
	goto st7;
tr129:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st114;
st114:
	if ( ++p == pe )
		goto _test_eof114;
goto case; case 114:
#line 3375 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr112;
	goto st7;
tr9:
#line 54 "sam_alignment.rl"
	{ rname_beg = p - line.ptr; }
	goto st115;
st115:
	if ( ++p == pe )
		goto _test_eof115;
goto case; case 115:
#line 3387 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr130;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st115;
	goto st5;
tr6:
#line 24 "sam_alignment.rl"
	{ int_value = 0; }
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st116;
st116:
	if ( ++p == pe )
		goto _test_eof116;
goto case; case 116:
#line 3403 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr133;
	goto st3;
tr133:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st117;
st117:
	if ( ++p == pe )
		goto _test_eof117;
goto case; case 117:
#line 3417 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr134;
	goto st3;
tr134:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st118;
st118:
	if ( ++p == pe )
		goto _test_eof118;
goto case; case 118:
#line 3431 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr135;
	goto st3;
tr135:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st119;
st119:
	if ( ++p == pe )
		goto _test_eof119;
goto case; case 119:
#line 3445 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr136;
	goto st3;
tr136:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st120;
st120:
	if ( ++p == pe )
		goto _test_eof120;
goto case; case 120:
#line 3459 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr137;
	goto st3;
tr137:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st121;
st121:
	if ( ++p == pe )
		goto _test_eof121;
goto case; case 121:
#line 3473 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr138;
	goto st3;
tr138:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st122;
st122:
	if ( ++p == pe )
		goto _test_eof122;
goto case; case 122:
#line 3487 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr139;
	goto st3;
tr139:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st123;
st123:
	if ( ++p == pe )
		goto _test_eof123;
goto case; case 123:
#line 3501 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr140;
	goto st3;
tr140:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st124;
st124:
	if ( ++p == pe )
		goto _test_eof124;
goto case; case 124:
#line 3515 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr141;
	goto st3;
tr141:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st125;
st125:
	if ( ++p == pe )
		goto _test_eof125;
goto case; case 125:
#line 3529 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr142;
	goto st3;
tr142:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st126;
st126:
	if ( ++p == pe )
		goto _test_eof126;
goto case; case 126:
#line 3543 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr143;
	goto st3;
tr143:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st127;
st127:
	if ( ++p == pe )
		goto _test_eof127;
goto case; case 127:
#line 3557 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr144;
	goto st3;
tr144:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st128;
st128:
	if ( ++p == pe )
		goto _test_eof128;
goto case; case 128:
#line 3571 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr145;
	goto st3;
tr145:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st129;
st129:
	if ( ++p == pe )
		goto _test_eof129;
goto case; case 129:
#line 3585 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr146;
	goto st3;
tr146:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st130;
st130:
	if ( ++p == pe )
		goto _test_eof130;
goto case; case 130:
#line 3599 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr147;
	goto st3;
tr147:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st131;
st131:
	if ( ++p == pe )
		goto _test_eof131;
goto case; case 131:
#line 3613 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr148;
	goto st3;
tr148:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st132;
st132:
	if ( ++p == pe )
		goto _test_eof132;
goto case; case 132:
#line 3627 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr149;
	goto st3;
tr149:
#line 25 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st133;
st133:
	if ( ++p == pe )
		goto _test_eof133;
goto case; case 133:
#line 3641 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr132;
	goto st3;
tr2:
#line 43 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
	goto st134;
st134:
	if ( ++p == pe )
		goto _test_eof134;
goto case; case 134:
#line 3653 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr150;
	if ( (*p) > 63u ) {
		if ( 65u <= (*p) && (*p) <= 126u )
			goto st134;
	} else if ( (*p) >= 33u )
		goto st134;
	goto st1;
		default: break;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof153: cs = 153; goto _test_eof; 
	_test_eof154: cs = 154; goto _test_eof; 
	_test_eof155: cs = 155; goto _test_eof; 
	_test_eof156: cs = 156; goto _test_eof; 
	_test_eof157: cs = 157; goto _test_eof; 
	_test_eof158: cs = 158; goto _test_eof; 
	_test_eof159: cs = 159; goto _test_eof; 
	_test_eof160: cs = 160; goto _test_eof; 
	_test_eof161: cs = 161; goto _test_eof; 
	_test_eof162: cs = 162; goto _test_eof; 
	_test_eof163: cs = 163; goto _test_eof; 
	_test_eof164: cs = 164; goto _test_eof; 
	_test_eof165: cs = 165; goto _test_eof; 
	_test_eof166: cs = 166; goto _test_eof; 
	_test_eof167: cs = 167; goto _test_eof; 
	_test_eof168: cs = 168; goto _test_eof; 
	_test_eof169: cs = 169; goto _test_eof; 
	_test_eof170: cs = 170; goto _test_eof; 
	_test_eof171: cs = 171; goto _test_eof; 
	_test_eof172: cs = 172; goto _test_eof; 
	_test_eof173: cs = 173; goto _test_eof; 
	_test_eof174: cs = 174; goto _test_eof; 
	_test_eof175: cs = 175; goto _test_eof; 
	_test_eof176: cs = 176; goto _test_eof; 
	_test_eof177: cs = 177; goto _test_eof; 
	_test_eof178: cs = 178; goto _test_eof; 
	_test_eof179: cs = 179; goto _test_eof; 
	_test_eof180: cs = 180; goto _test_eof; 
	_test_eof181: cs = 181; goto _test_eof; 
	_test_eof182: cs = 182; goto _test_eof; 
	_test_eof183: cs = 183; goto _test_eof; 
	_test_eof184: cs = 184; goto _test_eof; 
	_test_eof185: cs = 185; goto _test_eof; 
	_test_eof186: cs = 186; goto _test_eof; 
	_test_eof187: cs = 187; goto _test_eof; 
	_test_eof188: cs = 188; goto _test_eof; 
	_test_eof189: cs = 189; goto _test_eof; 
	_test_eof190: cs = 190; goto _test_eof; 
	_test_eof191: cs = 191; goto _test_eof; 
	_test_eof192: cs = 192; goto _test_eof; 
	_test_eof193: cs = 193; goto _test_eof; 
	_test_eof194: cs = 194; goto _test_eof; 
	_test_eof195: cs = 195; goto _test_eof; 
	_test_eof196: cs = 196; goto _test_eof; 
	_test_eof197: cs = 197; goto _test_eof; 
	_test_eof198: cs = 198; goto _test_eof; 
	_test_eof199: cs = 199; goto _test_eof; 
	_test_eof200: cs = 200; goto _test_eof; 
	_test_eof201: cs = 201; goto _test_eof; 
	_test_eof202: cs = 202; goto _test_eof; 
	_test_eof203: cs = 203; goto _test_eof; 
	_test_eof204: cs = 204; goto _test_eof; 
	_test_eof205: cs = 205; goto _test_eof; 
	_test_eof206: cs = 206; goto _test_eof; 
	_test_eof207: cs = 207; goto _test_eof; 
	_test_eof208: cs = 208; goto _test_eof; 
	_test_eof209: cs = 209; goto _test_eof; 
	_test_eof210: cs = 210; goto _test_eof; 
	_test_eof211: cs = 211; goto _test_eof; 
	_test_eof212: cs = 212; goto _test_eof; 
	_test_eof213: cs = 213; goto _test_eof; 
	_test_eof214: cs = 214; goto _test_eof; 
	_test_eof215: cs = 215; goto _test_eof; 
	_test_eof216: cs = 216; goto _test_eof; 
	_test_eof217: cs = 217; goto _test_eof; 
	_test_eof218: cs = 218; goto _test_eof; 
	_test_eof219: cs = 219; goto _test_eof; 
	_test_eof220: cs = 220; goto _test_eof; 
	_test_eof221: cs = 221; goto _test_eof; 
	_test_eof222: cs = 222; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 143: 
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 186: 
#line 157 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]); 
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 183: 
#line 161 "sam_alignment.rl"
	{
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]);
        current_tagvalue.setHexadecimalFlag();
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 204: 
	case 205: 
	case 206: 
	case 207: 
	case 208: 
	case 209: 
	case 210: 
	case 211: 
	case 212: 
	case 213: 
	case 214: 
	case 215: 
	case 216: 
	case 217: 
	case 218: 
	case 219: 
	case 220: 
	case 221: 
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 127 "sam_alignment.rl"
	{ 
        if (int_value < 0) {
            if (int_value >= byte.min) {
                current_tagvalue = Value(to!byte(int_value));
            } else if (int_value >= short.min) {
                current_tagvalue = Value(to!short(int_value));
            } else if (int_value >= int.min) {
                current_tagvalue = Value(to!int(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                current_tagvalue = Value(to!ubyte(int_value));
            } else if (int_value <= ushort.max) {
                current_tagvalue = Value(to!ushort(int_value));
            } else if (int_value <= uint.max) {
                current_tagvalue = Value(to!uint(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 191: 
	case 194: 
	case 195: 
	case 198: 
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 153 "sam_alignment.rl"
	{ 
        current_tagvalue = Value(float_value);
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 149: 
	case 150: 
	case 151: 
	case 152: 
	case 153: 
	case 154: 
	case 155: 
	case 156: 
	case 157: 
	case 158: 
	case 159: 
	case 160: 
	case 161: 
	case 162: 
	case 163: 
	case 164: 
	case 165: 
	case 166: 
#line 26 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 180 "sam_alignment.rl"
	{
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
	case 171: 
	case 174: 
	case 175: 
	case 178: 
#line 34 "sam_alignment.rl"
	{ 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 193 "sam_alignment.rl"
	{ 
        outbuffer.write(float_value); 
    }
#line 197 "sam_alignment.rl"
	{
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }
#line 225 "sam_alignment.rl"
	{ builder.put(current_tag, current_tagvalue); }
	break;
#line 4051 "sam_alignment.d"
		default: break;
	}
	}

	}

#line 321 "sam_alignment.rl"

    auto read = BamRead(line[read_name_beg .. read_name_end], 
                        sequence,
                        cigar.data,
                        builder.data);

    if (qual_ptr !is null && qual_index == sequence.length) {
        read.phred_base_quality = qual_ptr[0 .. sequence.length];
    }

    read.flag = flag;
    read.mapping_quality = mapping_quality;
    read.position = pos - 1; // we use 0-based coordinates, not 1-based
    read.template_length = template_length;
    read.next_pos = mate_pos - 1; // also 0-based
    read.ref_id = ref_id;
    read.next_ref_id = mate_ref_id;

    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto header = new SamHeader("@SQ\tSN:20\tLN:1234567");
    auto alignment = parseAlignmentLine(line, header);
    assert(alignment.read_name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60032);
    assert(alignment.mapping_quality == 25);
    assert(alignment.next_pos == 60032);
    assert(alignment.ref_id == 0);
    assert(alignment.next_ref_id == 0);
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
    assert(to!char(alignment["XT"]) == 'U');

    import std.stdio;
    import bio.bam.serialization.sam;
    import bio.bam.reference;

    ReferenceSequenceInfo info;
    info.name = "20";
    info.length = 1234567;

    auto invalid_cigar_string = "1\t100\t20\t50000\t30\tMZABC\t=\t50000\t0\tACGT\t####";
    alignment = parseAlignmentLine(invalid_cigar_string, header);
    assert(equal(alignment.sequence(), "ACGT"));

    auto invalid_tag_and_qual = "2\t100\t20\t5\t40\t27M30X5D\t=\t3\t10\tACT\t !\n\tX1:i:7\tX3:i:zzz\tX4:i:5";
    alignment = parseAlignmentLine(invalid_tag_and_qual, header);
    assert(alignment.phred_base_quality == [255, 255, 255]); // i.e. invalid
    assert(to!ubyte(alignment["X1"]) == 7);
    assert(alignment["X3"].is_nothing);
    assert(to!ubyte(alignment["X4"]) == 5);

}
