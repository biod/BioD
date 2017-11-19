# BioD [![Build Status](https://travis-ci.org/biod/BioD.svg?branch=master)](https://travis-ci.org/biod/BioD)

[BioD](https://github.com/biod/BioD) is a fast and memory efficient bioinformatics library written in the [D programming language](http://dlang.org).

BioD aims to:

* Provide a platform for writing high-performance bioinformatics applications in D. BioD achieves this by:
  - automatic parallelization of tasks where possible for example reading and writing BAM files
  - reducing the GC overhead by avoiding unnecessary memory allocations
* Offer support for manipulating common biological data formats
* Write clear documented and maintainable codebase

# Usage

BioD can be installed via the [D package manager](https://code.dlang.org/packages/biod).

See the [examples directory](https://github.com/biod/BioD/tree/master/examples)
for examples and usage.

BioD is also a crucial part of the [sambamba](https://github.com/biod/sambamba) tool.

# BioD contributors and support

See [contributors](https://github.com/biod/BioD/graphs/contributors). For support
contact

* [Artem Tarasov](https://github.com/lomereiter)
* [Pjotr Prins](https://github.com/pjotrp)


# License

BioD is licensed under the liberal MIT (expat) [license](./LICENSE).
