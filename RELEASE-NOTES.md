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
