# BioD [![Build Status](https://travis-ci.org/biod/BioD.svg?branch=master)](https://travis-ci.org/biod/BioD) [![DUB Package](https://img.shields.io/badge/dub-v0.1.0-red.svg)](https://code.dlang.org/packages/biod)

[BioD](https://github.com/biod/BioD) is a fast and memory efficient bioinformatics library written in the [D programming language](http://dlang.org).

BioD aims to:

* Provide a platform for writing high-performance bioinformatics applications in D. BioD achieves this by:
  - automatic parallelization of tasks where possible for example reading and writing BAM files
  - reducing the GC overhead by avoiding unnecessary memory allocations
* Offer support for manipulating common biological data formats

# Install

The current default is to provide the path to the checked out repo to the D-compiler. For example
in sambamba we use

    DFLAGS = -wi -I. -IBioD -g

# Usage

See the [examples directory](https://github.com/biod/BioD/tree/master/examples)
for examples and usage.

BioD is also a crucial part of the [sambamba](https://github.com/biod/sambamba) tool.

# Contributing

Simply clone the repository on github and put in a pull request.

# BioD contributors and support

See
[contributors](https://github.com/biod/BioD/graphs/contributors). For
support use the [issue tracker](https://github.com/biod/BioD/issues) or contact

* [Pjotr Prins](https://github.com/pjotrp)
* [Artem Tarasov](https://github.com/lomereiter)
* [George Githinji](https://github.com/biorelated)

# License

BioD is licensed under the liberal MIT (expat) [license](./LICENSE).
