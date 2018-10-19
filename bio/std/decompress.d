/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.decompress;

/**
   Streaming line reader which can be used for gzipped files. Note the
   current edition still uses the garbage collector.

   Conversion can happen between different encodings, provided the
   line terminator is ubyte = '\n'. GzipbyLine logic is modeled on
   ByLineImpl and readln function from std.stdio.
*/

import std.conv;
import std.exception;
import std.file;
import std.stdio;
import std.zlib: UnCompress;

struct GzipbyLine(Char) {

  File f;
  UnCompress decompress;
  Char[] line;
  ubyte[] uncompressed_buf;
  uint _bufsize;

  this(string gzipfn, uint bufsize=1024) {
    enforce(gzipfn.isFile);
    f = File(gzipfn,"r");
    decompress = new UnCompress();
    _bufsize = bufsize;
  }

  @disable this(this); // disable copy semantics;

  int opApply(scope int delegate(ubyte[]) dg) {
    foreach (ubyte[] buffer; f.byChunk(_bufsize))
    {
      auto uncompressedData = decompress.uncompress(buffer);
      dg(cast(ubyte[])uncompressedData);
    }
    return 0;
  }
}

unittest {
  writeln("Testing GzipbyLine");
  uint lines = 0;
  uint chars = 0;
  foreach(ubyte[] s; GzipbyLine!ubyte("test/data/BXD_geno.txt.gz",32)) {
    // test file contains 147 lines 3622 characters
    // writeln(s);
    chars += s.length;
  }
  assert(chars == 3622,"size " ~ to!string(chars));
  // assert(lines == 147,"lines " ~ to!string(lines));
}
