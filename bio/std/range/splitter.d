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
    return 0;
  }
}

unittest {
  auto s = cast(ubyte[])"hello 1 2 \t3  4 \n";
  assert(array(SimpleSplitConv!(ubyte[])(s)) == ["hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"  hello, 1 2 \t3  4 \n")) == ["","hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);
}
