/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.range.splitter;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.stdio;

import std.range.primitives;

immutable ubyte[] SPLIT_ON = [ 0x20, 0x09, 0x0A, ';', ',' ];

/**
   SimpleSplitConv takes a range R (typically a text line) and splits
   it/tokenizes it on a list of characters. Essentially fields/tokens
   are split by tabs, semi-colons or comma's and spaces. This compares
   to C's strtok(str, ", \t;").

   This routine happens often in bioinformatics and is a replacement
   for the much unsafer C strtok.  This edition should also handle
   UTF.

   The default is to split on space, newline, tab, semi-colon and
   comma.
*/

struct SimpleSplitConv(R)
  if (isInputRange!R)
{
  R list, split_on;

  this(R range, R splits_on = cast(R)SPLIT_ON) {
    list = range;
    split_on = splits_on;
  }

  int opApply(scope int delegate(R) dg) {
    size_t start = 0;
    bool in_whitespace = false;
    foreach(size_t pos, c; list) {
      if (canFind(split_on,c)) { // hit split char
        if (!in_whitespace) { // emit
          auto token = list[start..pos];
          dg(token);
        }
        start = pos+1;
        in_whitespace = true;
      } else {
        in_whitespace = false;
      }
    }
    if (!in_whitespace) { // emit final
      auto token = list[start..$];
      dg(token);
    }
    return 0;
  }
}

unittest {
  auto s = cast(ubyte[])"hello 1 2 \t3  4 \n";
  assert(array(SimpleSplitConv!(ubyte[])(s)) == ["hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"  hello, 1 2 \t3  4 \n")) == ["","hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"chr1:55365,55365,1")) == ["chr1:55365","55365","1"]);
}

R[] fast_splitter(R)(R range, R splits_on = cast(R)SPLIT_ON) {
  R[] tokens = new R[range.length]; // pre-allocate optimistially
  auto j = 0, prev_j = 0;
  bool in_whitespace = false;
  auto token_num = 0;
  for (; j<range.length ;) {
    bool found = false;
    auto check = range[j];
    foreach (c ; splits_on) {
      if (c==check) {
        found = true;
        break;
      }
    }
    if (found) { // hit split char
      if (!in_whitespace && j>0) {
        tokens[token_num] = range[prev_j..j];
        token_num++;
      }
      prev_j = j+1;
      in_whitespace = true;
    }
    else {
      in_whitespace = false;
    }
    j++;
  }
  if (!in_whitespace) { // emit final
    tokens[token_num] = range[prev_j..$];
    token_num++;
  }
  tokens.length = token_num;
  return tokens;
}

unittest {
  auto s = "hello 1 2 \t3  4 \n";
  writeln(fast_splitter(s).map!"to!string(a)");
  for (int x = 0; x < 4_000_000; x++) {
    assert(fast_splitter(s) == ["hello", "1", "2", "3", "4"]);
    // assert(array(FastSplitConv!(ubyte[])(cast(ubyte[])"  hello, 1 2 \t3  4 \n")) == ["","hello","1","2","3","4"]);
    // assert(array(FastSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);
    // assert(array(FastSplitConv!(ubyte[])(cast(ubyte[])"chr1:55365,55365,1")) == ["chr1:55365","55365","1"]);
  }
  /*
    real    0m1.731s
    user    0m1.732s
    sys     0m0.000s

    real    0m2.675s
    user    0m2.676s
    sys     0m0.000s

    real    0m3.733s
    user    0m3.736s
    sys     0m0.000s
  */
}
