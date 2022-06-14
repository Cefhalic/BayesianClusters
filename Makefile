HEADERS = $(sort $(wildcard include/*.hpp) )

# Files for library
LIBRARY_SOURCES = $(sort $(wildcard src/*.cpp) $(wildcard src/**/*.cpp) )
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${LIBRARY_SOURCES})
# LIBRARY_FILE = libClusterize.so

# Files for executables
EXECUTABLE_SOURCES = $(sort $(wildcard src/*.cxx) )
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/bin/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,%.exe,${EXECUTABLE_SOURCES})

DOXYGEN = doxygen/html/index.html
DOCUMENTATION = documentation/OptimizingTheMaths.pdf

DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath})))

.PHONY: all _all clean _cleanall

default: all

clean: _cleanall
_cleanall:
	rm -rf obj doxygen ${EXECUTABLES} ${DOCUMENTATION}

all: _all
build: _all
buildall: _all
# _all: ${LIBRARY_FILE} ${EXECUTABLES}
_all: ${DOCUMENTATION} ${DOXYGEN} ${EXECUTABLES}

FLAGS = -g -std=c++11 -march=native -O3 -lm `root-config --glibs --cflags --libs` -lMathMore -flto -MMD -MP -lboost_program_options

${DIRECTORIES}:
	mkdir -p $@

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	@echo -e "*************************************************************************************************************************************\n* Building Object Files\n*************************************************************************************************************************************"
	g++ -c ${FLAGS} -Iinclude -fPIC $< -o $@

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	g++ -c ${FLAGS} -Iinclude -fPIC $< -o $@

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

# ${LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES}
# 	g++ -shared $^ ${FLAGS} -o $@

# ${EXECUTABLES}: %.exe: obj/bin/%.o 
# 	g++ $< -L. -lClusterize ${FLAGS} -o $@


${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	@echo -e "*************************************************************************************************************************************\n* Building Executable\n*************************************************************************************************************************************"
	g++ $^ ${FLAGS} -o $@

${DOXYGEN}: ${HEADERS} ${LIBRARY_SOURCES} ${EXECUTABLE_SOURCES}
	@echo -e "*************************************************************************************************************************************\n* Generating Doxygen Documentation\n*************************************************************************************************************************************"
	doxygen Doxyfile

${DOCUMENTATION}:
	@echo -e "*************************************************************************************************************************************\n* Generating Maths Documentation\n*************************************************************************************************************************************"
	pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
	pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
