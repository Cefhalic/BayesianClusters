HEADERS = $(sort $(wildcard include/*.hpp) )

# Files for library
LIBRARY_SOURCES = $(sort $(wildcard src/Utilities/*.cpp) $(wildcard src/BayesianClustering/*.cpp) )
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${LIBRARY_SOURCES})

PYTHON_SOURCES = $(sort $(wildcard src/PythonBindings/*.cpp) )
PYTHON_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${PYTHON_SOURCES})

PYTHON_LIBRARY_FILE = python/BayesianClustering.so

# Files for executables
EXECUTABLE_SOURCES = $(sort $(wildcard src/*.cxx) )
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/bin/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,%.exe,${EXECUTABLE_SOURCES})

DOXYGEN = documentation/SoftwareManual.pdf
DOCUMENTATION = documentation/OptimizingTheMaths.pdf

DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${PYTHON_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath}))) extern

LIBPYTHON = $(shell ${CONDA_PREFIX}/bin/python -c "from sys import version_info; print( f'python{version_info[0]}.{version_info[1]}' )" )
LIBBOOSTPYTHON = $(shell ${CONDA_PREFIX}/bin/python -c "from sys import version_info; print( f'boost_python{version_info[0]}{version_info[1]}' )" )

.PHONY: clean all help cpp doxygen docs verbose

default: cpp
verbose: cpp

clean:
	rm -rf obj .doxygen ${EXECUTABLES} ${PYTHON_LIBRARY_FILE} ${DOCUMENTATION} ${DOXYGEN}

all : cpp doxygen docs 

help:
	@echo -e "Makefile for the Bayesian clustering project.\nUsage:"
	@echo "  - make                - Build code"
	@echo "  - make help           - Display this help message"
	@echo "  - make clean          - Tidy all build products"
	@echo "  - make all            - Build code, generate doxygen documentation and produce PDFs of latex sources"
	@echo "  - make cpp            - Build code"
	@echo "  - make verbose        - Build code, echoing the full command"
	@echo "  - make doxygen        - Generate doxygen documentation"
	@echo "  - make docs           - Produce PDFs of latex sources"
	@echo

cpp: ${EXECUTABLES} ${PYTHON_LIBRARY_FILE}
doxygen: ${DOXYGEN} 
docs: ${DOCUMENTATION}

# Adding the includes and libs for the locally built deps first
FLAGS = -L${CONDA_PREFIX}/lib -Iinclude -I${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/boost   \
        -lgsl -lgslcblas -lboost_program_options -lm -lpthread  \
        -g -std=c++14 -march=native -O3 -MMD -MP \
        \
        -Wno-deprecated-declarations
# Hide the annoying boost auto_ptr=>unique_ptr warning     

PYTHONFLAGS = -I${CONDA_PREFIX}/include/${LIBPYTHON} -l${LIBBOOSTPYTHON} -l${LIBPYTHON}

ifeq (verbose, $(filter verbose,$(MAKECMDGOALS)))

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/PythonBindings/%.o : src/PythonBindings/%.cpp | $$(dir obj/lib/PythonBindings/%.o)
	test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	g++ $< -o $@ -c ${PYTHONFLAGS} ${FLAGS} -fPIC

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	g++ $^ -o $@ ${FLAGS}

${PYTHON_LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES} ${PYTHON_OBJECT_FILES}
	g++ $^ -o $@ -shared ${PYTHONFLAGS} ${FLAGS}
else

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	@echo "Building Object Files | g++ -c ... $< -o $@"
	@g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	@echo "Building Object Files | g++ -c ... $< -o $@"
	@g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/PythonBindings/%.o : src/PythonBindings/%.cpp | $$(dir obj/lib/PythonBindings/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	@echo "Building Object Files | g++ -c ... $< -o $@"		
	@g++ $< -o $@ -c ${PYTHONFLAGS} ${FLAGS} -fPIC

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	@echo "Building Executable   | g++ ... -o $@"
	@g++ $^ -o $@ ${FLAGS}

${PYTHON_LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES} ${PYTHON_OBJECT_FILES}
	@echo "Building Executable   | g++ ... -o $@"
	@g++ $^ -o $@ -shared ${PYTHONFLAGS} ${FLAGS}
endif



${DIRECTORIES}:
	@echo "Making directory      | mkdir -p $@"
	@mkdir -p $@

${DOXYGEN}: ${HEADERS} ${LIBRARY_SOURCES} ${EXECUTABLE_SOURCES}
	@echo "Generating Doxygen Documentation: doxygen Doxyfile ---> $@"
	@doxygen utilities/Doxyfile
	@make -C .doxygen/latex
	@cp .doxygen/latex/refman.pdf $@

${DOCUMENTATION}:
	@echo "Generating Maths Documentation: pdflatex ... documentation/OptimizingTheMaths ---> documentation/OptimizingTheMaths.pdf"
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
