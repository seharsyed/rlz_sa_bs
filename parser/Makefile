CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

.PHONY: all clean

all: rlz_lp

rlz_lp: rlz_lp.cpp rlz_lp.hpp parser.hpp
	$(CXX) $(CXXFLAGS) rlz_lp.cpp -o rlz_lp

clean:
	rm -f rlz_lp
