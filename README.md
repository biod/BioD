BioD
====

Library goals are:

* Support for most common formats of biological data
* Fast and low-memory processing, that includes:
    - Automatic parallelization where possible (e.g. BAM reading/writing)
    - Avoiding unnecessary allocations to reduce GC overhead
    - Spending hours with profiler to optimize hot paths ;-)
* Clear, documented, and maintainable codebase
* Being suitable for one-off scripts, comparable in size with those written in Python/Ruby/Perl
* Providing a platform for writing high-performance bioinformatics applications in D
