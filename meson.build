project('BioD', 'd',
    meson_version : '>=0.48',
    license : 'MIT',
    version : '0.2.3',
    default_options : ['buildtype=debugoptimized']
)

project_soversion = '0'

src_dir = include_directories('.')
pkgc = import('pkgconfig')

extra_dflags = []
if get_option('optimize_strong')
    if meson.get_compiler('d').get_id() == 'gcc'
        message('No custom strong optimization flags are defined for GDC yet.')
    endif
    extra_dflags = ['-O3', '-release', '-enable-inlining', '-boundscheck=off']
endif

#
# Dependencies
#
zlib_dep = dependency('zlib')

#
# Sources
#
undead_src = [
    'contrib/undead/cstream.d',
    'contrib/undead/doformat.d',
    'contrib/undead/internal/file.d',
    'contrib/undead/stream.d',
    'contrib/undead/utf.d'
]

biod_src = [
    'bio/core/base.d',
    'bio/core/bgzf/block.d',
    'bio/core/bgzf/chunk.d',
    'bio/core/bgzf/compress.d',
    'bio/core/bgzf/constants.d',
    'bio/core/bgzf/inputstream.d',
    'bio/core/bgzf/outputstream.d',
    'bio/core/bgzf/virtualoffset.d',
    'bio/core/call.d',
    'bio/core/decompress.d',
    'bio/core/genotype.d',
    'bio/core/kmer.d',
    'bio/core/region.d',
    'bio/core/sequence.d',
    'bio/core/tinymap.d',
    'bio/core/utils/algo.d',
    'bio/core/utils/bylinefast.d',
    'bio/core/utils/exception.d',
    'bio/core/utils/format.d',
    'bio/core/utils/memoize.d',
    'bio/core/utils/outbuffer.d',
    'bio/core/utils/range.d',
    'bio/core/utils/roundbuf.d',
    'bio/core/utils/stream.d',
    'bio/core/utils/switchendianness.d',
    'bio/core/utils/tmpfile.d',
    'bio/core/utils/zlib.d',
    'bio/std/experimental/hts/bam/header.d',
    'bio/std/experimental/hts/bam/reader.d',
    'bio/std/experimental/hts/bam/writer.d',
    'bio/std/experimental/hts/bgzf.d',
    'bio/std/experimental/hts/bgzf_writer.d',
    'bio/std/experimental/hts/constants.d',
    'bio/std/experimental/hts/hashing.d',
    'bio/std/experimental/hts/logger.d',
    'bio/std/experimental/hts/pileup.d',
    'bio/std/experimental/hts/reads.d',
    'bio/std/experimental/hts/unpack.d',
    'bio/std/file/fai.d',
    'bio/std/file/fasta.d',
    'bio/std/file/fastq.d',
    'bio/std/genotype/maf.d',
    'bio/std/genotype/snp.d',
    'bio/std/hts/bam/abstractreader.d',
    'bio/std/hts/bam/bai/bin.d',
    'bio/std/hts/bam/baifile.d',
    'bio/std/hts/bam/bai/indexing.d',
    'bio/std/hts/bam/baseinfo.d',
    'bio/std/hts/bam/cigar.d',
    'bio/std/hts/bam/constants.d',
    'bio/std/hts/bam/md/core.d',
    'bio/std/hts/bam/md/operation.d',
    'bio/std/hts/bam/md/parse.d',
    'bio/std/hts/bam/md/reconstruct.d',
    'bio/std/hts/bam/multireader.d',
    'bio/std/hts/bam/pileup.d',
    'bio/std/hts/bam/randomaccessmanager.d',
    'bio/std/hts/bam/read.d',
    'bio/std/hts/bam/reader.d',
    'bio/std/hts/bam/readrange.d',
    'bio/std/hts/bam/reference.d',
    'bio/std/hts/bam/referenceinfo.d',
    'bio/std/hts/bam/region.d',
    'bio/std/hts/bam/splitter.d',
    'bio/std/hts/bam/tagvalue.d',
    'bio/std/hts/bam/validation/alignment.d',
    'bio/std/hts/bam/validation/samheader.d',
    'bio/std/hts/bam/writer.d',
    'bio/std/hts/iontorrent/flowcall.d',
    'bio/std/hts/iontorrent/flowindex.d',
    'bio/std/hts/sam/header.d',
    'bio/std/hts/sam/reader.d',
    'bio/std/hts/sam/utils/fastrecordparser.d',
    'bio/std/hts/sam/utils/recordparser.d',
    'bio/std/hts/snpcallers/maq.d',
    'bio/std/hts/snpcallers/simple.d',
    'bio/std/hts/thirdparty/msgpack.d',
    'bio/std/hts/utils/array.d',
    'bio/std/hts/utils/graph.d',
    'bio/std/hts/utils/samheadermerger.d',
    'bio/std/hts/utils/value.d',
    'bio/std/maf/block.d',
    'bio/std/maf/parser.d',
    'bio/std/maf/reader.d',
    'bio/std/range/splitter.d',
    'bio/std/sff/constants.d',
    'bio/std/sff/index.d',
    'bio/std/sff/read.d',
    'bio/std/sff/reader.d',
    'bio/std/sff/readrange.d',
    'bio/std/sff/utils/roundup.d',
    'bio/std/sff/writer.d'
]

tests_src = [
    'test/read_bam_file.d',
    'test/unittests.d'
]

#
# Includes
#
install_subdir('bio/', install_dir: 'include/d/bio/')
install_subdir('contrib/undead', install_dir: 'include/d/bio/contrib/')

#
# Library and pkg-config
#
biod_lib = both_libraries('biod',
        [undead_src, biod_src],
        dependencies: [zlib_dep],
        install: true,
        version: meson.project_version(),
        soversion: project_soversion,
        d_args: extra_dflags
)
pkgc.generate(name: 'biod',
              libraries: biod_lib,
              subdirs: 'd/bio/',
              version: meson.project_version(),
              description: 'D library for computational biology and bioinformatics'
)

#
# Tests
#
if get_option('tests')
    biod_test_exe = executable('biod_test',
        [undead_src,
         biod_src,
         tests_src],
        dependencies: [zlib_dep],
        d_unittest: true,
        d_args: extra_dflags
    )
    test('biod_tests',
         biod_test_exe,
         args: ['--DRT-gcopt=gc:precise disable:1 cleanup:none'],
         workdir: join_paths(meson.source_root(), 'test'),
         is_parallel: false
    )
endif
