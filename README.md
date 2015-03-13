BioD
====
BioD is a fast bioinformatics library written in D programming language.
The main aims of BioD are:

* Support for reading and writing common biological data formats.
* Fast and low-memory processing, that includes:
    - Automatic parallelization (e.g. BAM reading/writing)
    - Reduce GC overhead by avoiding unnecessary allocations
    - Spending hours with profiler to optimize hot paths ;-)
* Clear, documented, and maintainable codebase
* Support for one-off scripts, comparable in size with those written in Python/Ruby/Perl
* Provide a platform for writing high-performance bioinformatics applications in D
