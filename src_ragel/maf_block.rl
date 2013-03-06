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
module bio.maf.parser;
import std.conv, std.array;
import bio.maf.block;

%%{
    machine maf_block;

    # Common utilities for parsing integers and floats
    action update_sign { current_sign = fc == '-' ? -1 : 1; }
    action init_integer { int_value = 0; }
    action consume_next_digit { int_value *= 10; int_value += fc - '0'; }
    action take_sign_into_account { int_value *= current_sign; current_sign = 1; }

    sign = [\-+];

    uint = ([0-9]{1,18}) > init_integer $ consume_next_digit ;
    int = (sign >update_sign)? uint % take_sign_into_account ;

    action mark_float_start { float_beg = p - line.ptr; }
    action update_float_value { 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }

    float = ((sign? ((digit* '.'? digit+ ([eE] sign? digit+)?) | "inf") ) | "nan")
                > mark_float_start % update_float_value ;
    # --------------------------------------------------------------------------

    action set_score { block.score = float_value; }
    action set_pass { block.pass = int_value; }
    # Alignment block line
    score_vp = "score=" float % set_score;
    pass_vp = "pass=" uint % set_pass ;
    ab_name_value_pair = score_vp | pass_vp ;
    alignment_block_line = 'a' (space+ ab_name_value_pair)* ;

    # Common
    action src_begin { src_beg = p - line.ptr; }
    action set_src { sequence.source = line[src_beg .. p - line.ptr]; }
    action set_start { sequence.start = int_value; }
    action set_size { sequence.size = int_value; }
    action set_strand { sequence.strand = fc; }
    action set_src_size { sequence.source_size = int_value; }
    action add_sequence { sequences.put(sequence); sequence = MafSequence.init; }
    action check_sequence { assert(line[src_beg .. p - line.ptr] == sequences.data.back.source); }
    src = (^space)+ > src_begin % set_src ;
    start = uint % set_start ;
    size = uint % set_size ;
    strand = ('+' | '-') > set_strand ;
    srcSize = uint % set_src_size ; 

    # Sequence line
    action text_begin { text_beg = p - line.ptr; }
    action set_text { sequence.text = line[text_beg .. p - line.ptr]; }
    text = (^space)+ ;
    s_line = ('s'
              space+ src 
              space+ start 
              space+ size 
              space+ strand
              space+ srcSize
              space+ text > text_begin % set_text) % add_sequence ;

    # 'i' line
    action set_left_status { sequences.data.back.left_status = fc; }
    action set_left_count { sequences.data.back.left_count = int_value; }
    action set_right_status { sequences.data.back.right_status = fc; }
    action set_right_count { sequences.data.back.right_count = int_value; }
    i_status = [CINnMT] ;
    leftStatus = i_status ;
    leftCount = uint ;
    rightStatus = i_status ;
    rightCount = uint ;
    i_line = 'i'
             space+ (src > src_begin % check_sequence)
             space+ leftStatus > set_left_status
             space+ leftCount % set_left_count
             space+ rightStatus > set_right_status
             space+ rightCount % set_right_count ;

    # 'e' line
    action set_empty_status { sequence.empty_status = *p; }
    e_status = [CIMn] ;
    e_line = ('e'
              space+ src
              space+ start
              space+ size
              space+ strand
              space+ srcSize
              space+ (e_status > set_empty_status)) % add_sequence ;

    # 'q' line
    action qual_begin { qual_beg = p - line.ptr; }
    action set_qual { sequences.data.back.quality = line[qual_beg .. p - line.ptr]; }
    q_value = (digit | 'F' | '-')+ ;
    q_line = 'q'
             space+ (src > src_begin % check_sequence)
             space+ (q_value > qual_begin % set_qual);

    newline = "\n" | "\r\n" ;
    block := alignment_block_line space* 
             (newline ((s_line | i_line | e_line | q_line) space*))+ ;

    write data; 
}%%

MafBlock parseMafBlock(string line) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    int current_sign;
    int int_value;
    double float_value;
    size_t float_beg;

    MafBlock block;
    MafSequence sequence;
    auto sequences = Appender!(MafSequence[])();

    size_t src_beg;
    size_t text_beg;
    size_t qual_beg;

    %%write init;
    %%write exec;

    block.sequences = sequences.data;
    return block;
}
