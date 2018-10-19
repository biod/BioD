# Simple Makefile

D_COMPILER=ldc2
DFLAGS = -wi -g -relocation-model=pic -unittest -main

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

SRC         = $(wildcard bio/std/*.d) bio2/bgzf.d bio/core/bgzf/constants.d
OBJ         = $(SRC:.d=.o)
BIN         = bin/biod_tests

debug:          DFLAGS += -O0 -d-debug -link-debuglib
release static: DFLAGS += -O3 -release -enable-inlining -boundscheck=off
static:         DFLAGS += -static -L-Bstatic

all: debug

default: all

default debug release static: $(BIN)

%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

$(BIN): $(OBJ)
	$(info linking...)
	$(D_COMPILER) $(DFLAGS) $(OBJ) -of=$(BIN)

check: $(BIN)
	$(BIN)
clean:
	rm -vf $(OBJ)
	rm -v $(BIN)
