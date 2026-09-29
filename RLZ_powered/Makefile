
ifndef LARGE_BLOCK_SIZE
LARGE_BLOCK_SIZE = 16384
endif

ifndef SMALL_BLOCK_SIZE
SMALL_BLOCK_SIZE = 256
endif

CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native -fopenmp\
         -DLARGE_BLOCK_SIZE=$(LARGE_BLOCK_SIZE) \
         -DSMALL_BLOCK_SIZE=$(SMALL_BLOCK_SIZE) \
         -DRUN_COUNT=$(RUN_COUNT) \
         -DCACHE_LINE=$(CL)

# -------- user parameters --------
INPUT ?=
OUTPUT ?=
FLAGS ?=

# -------- headers --------
HEADERS = include/reader.hpp include/block_rlbwt.hpp \
          include/byte_alphabet.hpp include/super_block.hpp \
          include/types.hpp include/two_byte_block.hpp include/custom_alphabet.hpp \
          include/alphabet.hpp include/non_rle_block.hpp

.PHONY: clean update_git debug all

.DEFAULT: all

%/%.hpp:

# -------- default target --------
all: gpp make_alphabet_header
gpp: make_bwt rlz_parser decompress compare_files

# -------- build binaries --------
make_bwt: make_bwt.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -o make_bwt make_bwt.cpp

make_alphabet_header: make_alphabet_header.cpp include/reader.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -o make_alphabet_header make_alphabet_header.cpp

rlz_parser: rlz_parser.cpp $(HEADERS)
	g++ $(CFLAGS) -DNDEBUG -Ofast -o rlz_parser rlz_parser.cpp

decompress: decompress.cpp
	g++ $(CFLAGS) -Ofast -o decompress decompress.cpp
	
compare_files: compare_files.cpp $(HEADERS)
	g++ $(CFLAGS) -Ofast -o compare_files compare_files.cpp


# -------- debug build --------
debug:
	g++ $(CFLAGS) -DDEBUG -g -o make_bwt make_bwt.cpp
	g++ $(CFLAGS) -DDEBUG -g -o rlz_parser rlz_parser.cpp

# -------- clean --------
clean:
	rm -f make_bwt rlz_parser make_alphabet_header decompress compare_files
