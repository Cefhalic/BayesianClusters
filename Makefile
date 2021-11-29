

.PHONY: all _all # clean _cleanall

default: all

# clean: _cleanall
# _cleanall:
# 	rm -rf obj lib bin sbin elements

all: _all
build: _all
buildall: _all
_all: test.exe


test.exe : src/test.cxx
	g++ -march=native -std=c++11 -Wall -O3 $^ -lm `root-config --glibs --cflags --libs` -o $@
