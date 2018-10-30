project('BioD', 'd',
    meson_version : '>=0.45',
    license : 'MIT',
    version : '0.2.2'
)

project_soversion = '0'

src_dir = include_directories('.')
pkgc = import('pkgconfig')

#
# Dependencies
#
# undead_dep = dependency('undead', method: 'dub')
undead_path = '/home/pjotr/.dub/packages/undead-1.0.9/undead/src'
undead_dir = include_directories('/home/pjotr/.dub/packages/undead-1.0.9/undead/src')

zlib_dep = dependency('zlib')

# thread_dep = dependency('threads')

#
# Sources
#
biod_src = [
    'bio/bam/abstractreader.d',
    'bio/bam/bai/bin.d',
    'bio/bam/baifile.d',
    'bio/bam/bai/indexing.d',
    'bio/bam/baseinfo.d',
    'bio/bam/cigar.d',
    'bio/bam/constants.d',
    'bio/bam/iontorrent/flowcall.d',
    'bio/bam/iontorrent/flowindex.d',
    'bio/bam/md/core.d',
    'bio/bam/md/operation.d',
    'bio/bam/md/parse.d',
    'bio/bam/md/reconstruct.d',
    'bio/bam/multireader.d',
    'bio/bam/pileup.d',
    'bio/bam/randomaccessmanager.d',
    'bio/bam/read.d',
    'bio/bam/reader.d',
    'bio/bam/readrange.d',
    'bio/bam/reference.d',
    'bio/bam/referenceinfo.d',
    'bio/bam/region.d',
    'bio/bam/snpcallers/maq.d',
    'bio/bam/snpcallers/simple.d',
    'bio/bam/splitter.d',
    'bio/bam/tagvalue.d',
    'bio/bam/thirdparty/msgpack.d',
    'bio/bam/utils/array.d',
    'bio/bam/utils/graph.d',
    'bio/bam/utils/samheadermerger.d',
    'bio/bam/utils/value.d',
    'bio/bam/validation/alignment.d',
    'bio/bam/validation/samheader.d',
    'bio/bam/writer.d',
    'bio/core/base.d',
    'bio/core/bgzf/block.d',
    'bio/core/bgzf/chunk.d',
    'bio/core/bgzf/compress.d',
    'bio/core/bgzf/constants.d',
    'bio/core/bgzf/inputstream.d',
    'bio/core/bgzf/outputstream.d',
    'bio/core/bgzf/virtualoffset.d',
    'bio/core/call.d',
    'bio/core/fasta.d',
    'bio/core/fastq.d',
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
    'bio/maf/block.d',
    'bio/maf/parser.d',
    'bio/maf/reader.d',
    'bio/sam/header.d',
    'bio/sam/reader.d',
    'bio/sam/utils/fastrecordparser.d',
    'bio/sam/utils/recordparser.d',
    'bio/sff/constants.d',
    'bio/sff/index.d',
    'bio/sff/read.d',
    'bio/sff/reader.d',
    'bio/sff/readrange.d',
    'bio/sff/utils/roundup.d',
    'bio/sff/writer.d',
    'bio/std/decompress.d',
    'bio/std/range/splitter.d'
]

bio2_src = [
    'bio2/bam/header.d',
    'bio2/bam/reader.d',
    'bio2/bam/writer.d',
    'bio2/bgzf.d',
    'bio2/bgzf_writer.d',
    'bio2/constants.d',
    'bio2/hashing.d',
    'bio2/logger.d',
    'bio2/pileup.d',
    'bio2/reads.d',
    'bio2/unpack.d'
]

undead_src = [
    undead_path+'/undead/bitarray.d',
    undead_path+'/undead/cstream.d',
    undead_path+'/undead/datebase.d',
    undead_path+'/undead/date.d',
    undead_path+'/undead/dateparse.d',
    undead_path+'/undead/doformat.d',
    undead_path+'/undead/metastrings.d',
    undead_path+'/undead/regexp.d',
    undead_path+'/undead/socketstream.d',
    undead_path+'/undead/stream.d',
    undead_path+'/undead/string.d',
    undead_path+'/undead/utf.d',
    undead_path+'/undead/internal/file.d'
]

#
# Includes
#
install_subdir('bio/', install_dir: 'include/d/bio/')
install_subdir('bio2/', install_dir: 'include/d/bio/')

add_global_arguments('-unittest', language : 'D')

#
# Library and pkg-config
#
biod_lib = library('biod',
        [biod_src, bio2_src, undead_src],
        include_directories : undead_dir,
        dependencies: [zlib_dep],
        install: true,
        d_unittest: true,
        link_args: ['-unittest'],
        version: meson.project_version(),
        soversion: project_soversion
)
pkgc.generate(name: 'biod',
              libraries: biod_lib,
              subdirs: 'd/bio/',
              version: meson.project_version(),
              description: 'Bioinformatics library in D (utils for working with SAM, BAM, SFF formats).'
)

#
# Tests
#

# if meson.get_compiler('d').get_id() == 'llvm'
#   extra_args = ['-unittest', '-link-defaultlib-shared']
# endif

biod_test_exe = executable('biod_test',
    ['test/unittests.d',
     biod_src,
     bio2_src, undead_src],
    include_directories : undead_dir,
    dependencies: [zlib_dep],
    d_unittest: true,
    c_args: ['-unittest'],
    link_args: ['-unittest']
)

test('biod_tests', biod_test_exe)
