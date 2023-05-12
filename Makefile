HEADERS = $(sort $(wildcard include/*.hpp) )

# Files for library
LIBRARY_SOURCES = $(sort $(wildcard src/Utilities/*.cpp) $(wildcard src/BayesianClustering/*.cpp) )
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${LIBRARY_SOURCES})

LIBRARY_FILE = libBayesianClusteringCore.so

PYTHON_SOURCES = $(sort $(wildcard src/PythonBindings/*.cpp) )
PYTHON_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${PYTHON_SOURCES})

PYTHON_LIBRARY_FILE = python/BayesianClustering.so

# Files for executables
EXECUTABLE_SOURCES = $(sort $(wildcard src/*.cxx) )
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/bin/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,%.exe,${EXECUTABLE_SOURCES})

DOXYGEN = documentation/SoftwareManual.pdf
DOCUMENTATION = documentation/OptimizingTheMaths.pdf

DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${PYTHON_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath})))

LIBPYTHON = $(shell ${CONDA_PREFIX}/bin/python -c "from sys import version_info; print( f'python{version_info[0]}.{version_info[1]}' )" )
LIBBOOSTPYTHON = $(shell ${CONDA_PREFIX}/bin/python -c "from sys import version_info; print( f'boost_python{version_info[0]}{version_info[1]}' )" )

FLAGS = -L${CONDA_PREFIX}/lib -Iinclude -I${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/boost   \
        -lgsl -lgslcblas -lboost_program_options -lm -lpthread  \
        -g -std=c++14 -march=native -O3 -MMD -MP -fPIC
      
PYTHONFLAGS = -I${CONDA_PREFIX}/include/${LIBPYTHON} -l${LIBBOOSTPYTHON} -l${LIBPYTHON} \
              -Wno-deprecated-declarations # Hide the annoying boost auto_ptr=>unique_ptr warning     


ifeq (verbose, $(filter verbose,$(MAKECMDGOALS)))
define switch_verbose 
	${2}
endef
else
define switch_verbose 
	@echo ${1}; ${2}
endef
endif


.PHONY: clean all help cpp doxygen docs verbose

default: cpp
verbose: cpp

clean:
	rm -rf obj .doxygen ${LIBRARY_FILE} ${EXECUTABLES} ${PYTHON_LIBRARY_FILE}

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

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	$(call switch_verbose, "Building Object Files | g++ -c ... $< -o $@" , g++ $< -o $@ -c ${FLAGS} )

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	$(call switch_verbose, "Building Object Files | g++ -c ... $< -o $@" , g++ $< -o $@ -c ${FLAGS} )

.SECONDEXPANSION:
obj/lib/PythonBindings/%.o : src/PythonBindings/%.cpp | $$(dir obj/lib/PythonBindings/%.o)
	@test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	$(call switch_verbose, "Building Object Files | g++ -c ... $< -o $@" , g++ $< -o $@ -c ${PYTHONFLAGS} ${FLAGS} )

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(PYTHON_OBJECT_FILES:.o=.d)	
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES}
	$(call switch_verbose, "Building Library      | g++ ... -o $@" , g++ $^ -o $@ -shared ${FLAGS} )

${PYTHON_LIBRARY_FILE}: ${LIBRARY_FILE} ${PYTHON_OBJECT_FILES}
	$(call switch_verbose, "Building Library      | g++ ... -o $@" , g++ $^ -o $@ -shared -L. -lBayesianClusteringCore ${PYTHONFLAGS} ${FLAGS} )

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_FILE}
	$(call switch_verbose, "Building Executable   | g++ ... -o $@" , g++ $^ -o $@ -L. -lBayesianClusteringCore ${FLAGS} )

${DIRECTORIES}:
	@echo "Making directory      | mkdir -p $@"
	@mkdir -p $@

${DOXYGEN}: ${HEADERS} ${LIBRARY_SOURCES} ${EXECUTABLE_SOURCES}
	@echo "Generating Doxygen Documentation: doxygen utilities/Doxyfile ---> $@"
	@doxygen utilities/Doxyfile
	@make -sC .doxygen/latex > /dev/null 2>&1
	@cp .doxygen/latex/refman.pdf $@

${DOCUMENTATION}:
	@echo "Generating Maths Documentation: pdflatex ... documentation/OptimizingTheMaths ---> documentation/OptimizingTheMaths.pdf"
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
