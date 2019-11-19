# Simple Makefile
#
#   make sharedlibrary  : make shared library

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

SRC         = $(wildcard contrib/undead/*.d) contrib/undead/*/*.d $(wildcard bio/*.d bio/*/*.d bio/*/*/*.d bio/*/*/*/*.d bio/*/*/*/*/*.d bio/*/*/*/*/*/*.d) test/unittests.d
OBJ         = $(SRC:.d=.o)
BIN         = bin/biod_tests
sharedlibrary:  BIN = libbiod.so

debug check:    DFLAGS += -O0 -d-debug -unittest -link-debuglib
release static: DFLAGS += -O3 -release -enable-inlining -boundscheck=off
static:         DFLAGS += -static -L-Bstatic
sharedlibrary:  DFLAGS += -shared

all: debug

default: all

default debug release static sharedlibrary: $(BIN)

%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

$(BIN): $(OBJ)
	$(info linking...)
	$(D_COMPILER) -main $(DFLAGS) $(OBJ) -of=$(BIN)

check: $(BIN)
	$(info running tests...)
	$(BIN) "--DRT-gcopt=gc:precise disable:1 cleanup:none"

clean:
	rm -vf $(OBJ)
	rm -v $(BIN)
        # find -name '*.o' -exec rm \{\} \;
