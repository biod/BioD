# BioD [![Build Status](https://travis-ci.org/biod/BioD.svg?branch=master)](https://travis-ci.org/biod/BioD) [![DUB Package](https://img.shields.io/badge/dub-v0.1.0-red.svg)](https://code.dlang.org/packages/biod)

[BioD](https://github.com/biod/BioD) is a fast and memory efficient bioinformatics library written in the [D programming language](http://www.dlang.org)
whose aim is to:

* Provide a platform for developing high-performance computational biology applications using the [D programming language](http://www.dlang.org) through
  - Automatic parallelization of tasks where possible
  - Avoiding unnecessary memory allocations

## Why BioD?

BioD leverages on [D programming language](http://www.dlang.org)
features to develop high performance bioinformatics tools
(e.g. [sambamba](https://github.com/biod/sambamba)). The D programming
language is both a low and high-level hybrid object orientated and
functional (OOP/FP) programming language with templating/generic
features are far easier than that of C++.

## D programming language resources

* [Programming in D](http://ddili.org/ders/d.en/index.html) is online by Ali Çehreli.
* [The D Programming Language](https://www.amazon.com/D-Programming-Language-Andrei-Alexandrescu/dp/0321635361) by Andrei Alexandrecu (great book, slightly out of date)
* [The D Cookbook](https://www.amazon.com/D-Cookbook-Adam-D-Ruppe/dp/1783287217) by Adam D. Ruppe

## Current development

Our aim is to provide a set of D modules to manipulate and work with
biological datasets.  BioD provides modules for manipulating high
throughput data formats by provifing fast and easy to use native BAM
file reader and writer with ability to iterate a BAM file a read at a
time,a nucleotide at a time (pileup) or via a sliding window.

Note the current Bamreader bails out on recent versions of the LDC
compiler. See also https://github.com/biod/BioD/issues/53

## Install

The current default is to provide the path to the checked out repo to
the D-compiler. For example,

    DFLAGS = -wi -I. -IBioD -g

## Build environment

After installing ldc and dub

    dub
    dub test

On a recent Debian (>201911) you can install ldc and compile BioD with

    make
    make check

It is possible to create a recent build container with the
[GNU guix](https://www.gnu.org/software/guix/) transactional package
manager

    guix environment -C guix --ad-hoc ldc dub zlib gdb binutils-gold vim --network

after getting dropped in the container simply run dub.

If you want to use the make file instead (not requiring the network) use

    guix environment -C guix --ad-hoc ldc zlib gdb make binutils-gold vim --no-grafts
    make -j 4
    make check

## Debugging

When using gdb, switch off these handlers

    handle SIGUSR1 SIGUSR2 nostop noprint

It can be passed in from the command line

    gdb -ex 'handle SIGUSR1 SIGUSR2 nostop noprint' --args ./biod-test-library

## Usage

See the [examples directory](examples)
for examples and usage.

## Mailing list

[The BioD mailing list](https://groups.google.com/forum/#!forum/dlang_biod)

## Contributing

Simply clone the repository on github and put in a pull request.

## BioD contributors and support

See
[contributors](https://github.com/biod/BioD/graphs/contributors). For
support use the [issue tracker](https://github.com/biod/BioD/issues) or contact

* [Pjotr Prins](https://github.com/pjotrp)
* [George Githinji](https://github.com/George-Githinji)
* [Prasun Anand](https://github.com/prasunanand)

## License

BioD is free software and licensed under the MIT (expat) [license](./LICENSE).

BioD includes some files from the
[undeaD project](https://github.com/dlang/undeaD) in ./contrib which
are published under a Boost license. This code should be phased out in time.
