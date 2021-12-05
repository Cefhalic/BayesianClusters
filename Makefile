# Files for library
LIBRARY_SOURCES = $(sort $(wildcard src/*.cpp) $(wildcard src/**/*.cpp) )
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${LIBRARY_SOURCES})
LIBRARY_FILE = libClusterize.so

# Files for executables
EXECUTABLE_SOURCES = $(sort $(wildcard src/*.cxx) )
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/bin/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,%.exe,${EXECUTABLE_SOURCES})


DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath})))

.PHONY: all _all # clean _cleanall

default: all

# clean: _cleanall
# _cleanall:
# 	rm -rf obj lib bin sbin elements

all: _all
build: _all
buildall: _all
_all: ${LIBRARY_FILE} ${EXECUTABLES}


# test.exe : src/test.cxx
# 	g++ -march=native -std=c++11 -Wall -O3 $^ -lm `root-config --glibs --cflags --libs` -lMathMore -o $@

FLAGS = -g -std=c++11 -march=native -O3 -lm `root-config --glibs --cflags --libs` -lMathMore


${DIRECTORIES}:
	mkdir -p $@


.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	g++ -c -march=native -std=c++11 -Wall -O3 -lm `root-config --glibs --cflags --libs` -lMathMore -Iinclude -fPIC $< -o $@

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	g++ -c -march=native -std=c++11 -Wall -O3 -lm `root-config --glibs --cflags --libs` -lMathMore -Iinclude -fPIC $< -o $@


-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)


${LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES}
	g++ -shared $^ ${FLAGS} -o $@

${EXECUTABLES}: %.exe: obj/bin/%.o 
	g++ $< -L. -lClusterize ${FLAGS} -o $@
