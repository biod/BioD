project('BioD', 'd')

project_version   = '0.1.0'
project_soversion = '0'

src_dir = include_directories('.')
pkgc = import('pkgconfig')

#
# Dependencies
#
undead_dep = dependency('undead', version : '>=1.0.6')
zlib_dep = dependency('zlib')

#
# Sources
#
biod_src = [
    'bio/bam/utils/value.d',
    'bio/bam/utils/samheadermerger.d',
    'bio/bam/utils/array.d',
    'bio/bam/utils/graph.d',
    'bio/bam/pileup.d',
    'bio/bam/readrange.d',
    'bio/bam/reader.d',
    'bio/bam/thirdparty/msgpack.d',
    'bio/bam/baifile.d',
    'bio/bam/iontorrent/flowcall.d',
    'bio/bam/iontorrent/flowindex.d',
    'bio/bam/splitter.d',
    'bio/bam/abstractreader.d',
    'bio/bam/bai/indexing.d',
    'bio/bam/bai/bin.d',
    'bio/bam/multireader.d',
    'bio/bam/reference.d',
    'bio/bam/md/operation.d',
    'bio/bam/md/core.d',
    'bio/bam/md/parse.d',
    'bio/bam/md/reconstruct.d',
    'bio/bam/tagvalue.d',
    'bio/bam/validation/alignment.d',
    'bio/bam/validation/samheader.d',
    'bio/bam/region.d',
    'bio/bam/randomaccessmanager.d',
    'bio/bam/referenceinfo.d',
    'bio/bam/constants.d',
    'bio/bam/snpcallers/maq.d',
    'bio/bam/snpcallers/simple.d',
    'bio/bam/read.d',
    'bio/bam/writer.d',
    'bio/bam/baseinfo.d',
    'bio/maf/reader.d',
    'bio/maf/parser.d',
    'bio/maf/block.d',
    'bio/core/utils/memoize.d',
    'bio/core/utils/algo.d',
    'bio/core/utils/stream.d',
    'bio/core/utils/zlib.d',
    'bio/core/utils/tmpfile.d',
    'bio/core/utils/roundbuf.d',
    'bio/core/utils/format.d',
    'bio/core/utils/bylinefast.d',
    'bio/core/utils/range.d',
    'bio/core/utils/switchendianness.d',
    'bio/core/utils/outbuffer.d',
    'bio/core/genotype.d',
    'bio/core/fasta.d',
    'bio/core/base.d',
    'bio/core/kmer.d',
    'bio/core/region.d',
    'bio/core/tinymap.d',
    'bio/core/sequence.d',
    'bio/core/bgzf/outputstream.d',
    'bio/core/bgzf/chunk.d',
    'bio/core/bgzf/inputstream.d',
    'bio/core/bgzf/constants.d',
    'bio/core/bgzf/block.d',
    'bio/core/bgzf/virtualoffset.d',
    'bio/core/bgzf/compress.d',
    'bio/core/call.d',
    'bio/sff/utils/roundup.d',
    'bio/sff/index.d',
    'bio/sff/readrange.d',
    'bio/sff/reader.d',
    'bio/sff/constants.d',
    'bio/sff/read.d',
    'bio/sff/writer.d',
    'bio/sam/utils/fastrecordparser.d',
    'bio/sam/utils/recordparser.d',
    'bio/sam/reader.d',
    'bio/sam/header.d'
]

#
# Includes
#
install_subdir('bio/', install_dir: 'include/d/')

#
# Library and pkg-config
#
biod_lib = library('biod',
        [biod_src],
        include_directories: [src_dir],
        dependencies: [undead_dep, zlib_dep],
        install: true,
        version: project_version,
        soversion: project_soversion
)
pkgc.generate(name: 'biod',
              libraries: biod_lib,
              subdirs: 'd/',
              version: project_version,
              description: 'Bioinformatics library in D (utils for working with SAM, BAM, SFF formats).'
)

#
# Tests
#
biod_test_exe = executable('biod_test',
    ['test/unittests.d',
    biod_src],
    include_directories: [src_dir],
    dependencies: [undead_dep, zlib_dep],
    d_args: meson.get_compiler('d').unittest_args()
)
test('biod_tests', biod_test_exe)
