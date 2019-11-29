# Simple Makefile
#
#   make shared  : make shared lib
#   make lib     : make static lib (nyi)
#   make check
#
# You can also use 'dub' and 'dub test' instead

D_COMPILER=ldc2
DFLAGS = -wi -g -relocation-model=pic -Icontrib/undead -L-lz

ifndef GUIX
  ifdef GUIX_ENVIRONMENT
    GUIX=$(GUIX_ENVIRONMENT)
  endif
endif
ifdef GUIX
  LIBRARY_PATH=$(GUIX)/lib
endif

DLIBS       = $(LIBRARY_PATH)/libphobos2-ldc.a $(LIBRARY_PATH)/libdruntime-ldc.a
DLIBS_DEBUG = $(LIBRARY_PATH)/libphobos2-ldc-debug.a $(LIBRARY_PATH)/libdruntime-ldc-debug.a

SRC         = $(wildcard contrib/undead/*.d) contrib/undead/*/*.d $(wildcard bio/*.d bio/*/*.d bio/*/*/*.d bio/*/*/*/*.d bio/*/*/*/*/*.d bio/*/*/*/*/*/*.d) test/unittests.d test/read_bam_file.d

OBJ         = $(SRC:.d=.o)
BIN         = bin/biod_tests
shared:     LIB = libbiod.so
lib:        LIB = libbiod

# debug check:    DFLAGS += -O0 -d-debug -unittest -link-debuglib
check:          DFLAGS += -O0 -d-debug -link-debuglib -unittest
release static: DFLAGS += -O3 -release -enable-inlining -boundscheck=off
static:         DFLAGS += -static -L-Bstatic
shared:         DFLAGS += -shared
lib:            DFLAGS += -lib

all: debug

default: all

default debug release static: $(BIN)
shared lib: $(LIB)

%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

$(LIB): $(OBJ)
	$(info linking lib...)
	$(D_COMPILER) $(DFLAGS) $(OBJ) -of=$(LIB)

$(BIN): $(OBJ)
	$(info linking...)
	$(D_COMPILER) $(DFLAGS) $(OBJ) -of=$(BIN)

check: $(BIN)
	$(info Make check running tests...)
	$(BIN)

# $(BIN) "--DRT-gcopt=gc:precise disable:1 cleanup:none"

clean:
	rm -vf $(OBJ)
	rm -v $(BIN)
        # find -name '*.o' -exec rm \{\} \;
