## ChangeLog v0.2.4 (20191128)

+ Fixed major bug causing https://github.com/biod/sambamba/issues/393
+ Fixed dub and make files
+ Dub test still fails because of GC
+ Debian package on the way (thanks https://github.com/atille https://github.com/biod/BioD/issues/50)

## ChangeLog v0.2.3 (20191119)

+ Compiles and tests pass on Debian with dub and ldc 1.17.0

## ChangeLog v0.2.2 (20190316)

+ Restored make so we can compile without dub
+ Fasta and fasta indexing (.fai) support added (thanks Emilio Palumbo https://github.com/emi80)
+ Mate pair comparison and HI tag support added for BAM (thanks https://github.com/emi80)
+ Added Picard-style comparison for BAM (thanks https://github.com/TimurIs)
+ Added fast whitespace line splitter/tokenizer, a Phobos-style version and a faster C-style version (thanks https://github.com/pjotrp)
+ Added multi-allelic frequencies (MAF) support (thanks https://github.com/pjotrp)
+ Name spaces and directories reorganised (thanks George Gethinji https://github.com/george-githinji)
+ Pulled in D's undead repo (dropped dependency) and minimalised it to actual used files (@pjotrp)

## ChangeLog v0.2.1 (20181004)

+ Fix bunch of deprecation warnings
+ Fixes sort removing tags from @RG header https://github.com/biod/sambamba/issues/356
+ Fixes https://github.com/biod/BioD/issues/37

## ChangeLog v0.2.0 (20180915)

64-bit compilation should be fine on ldc 1.10

+ Added FastQ parser (thanks George Githinji @george-githinji)
+ Added rewritten Bamreader f4a6c1c55c8903e948b793adbb150704d1e267b2 (thanks Pjotr Prins @pjotrp)
+ Improved big-endian support and other bug fixes (thanks Artem Tarasov @lomereiter)
+ Include SAM/BAM AH tag (thanks Indraniel Das @indraniel)
+ D undeaD compilation fixes (thanks John Colvin @John-Colvin)
+ Meson build definition added (thanks Matthias Klumpp @ximion)
+ Update outputstream.d (thanks Brett T. Hannigan @godotgildor)
+ Moved Cigar into its own module 7994406592f277bc6950739aaf75c5f948cb7928
+ Added 'asserte' method which throws an exception on assert (failure)
+ Bug fixes:
  * #31 Bug in bio/core/utils/outbuffer.d: _heap.length is not correctly set affects reading long SAM records
  * https://github.com/biod/sambamba/issues/210 correct pseudo-bin calculation (thanks @jfarek)
  * JSON header output fixes #331
