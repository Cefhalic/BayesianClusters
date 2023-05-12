HEADERS = $(sort $(wildcard include/*.hpp) )

# Files for library
LIBRARY_SOURCES = $(sort $(wildcard src/*.cpp) $(wildcard src/**/*.cpp) )
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/lib/%.o,${LIBRARY_SOURCES})

PYTHON_LIBRARY_FILE = python/BayesianClustering.so

# Files for executables
EXECUTABLE_SOURCES = $(sort $(wildcard src/*.cxx) )
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/bin/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,%.exe,${EXECUTABLE_SOURCES})

DOXYGEN = documentation/SoftwareManual.pdf
DOCUMENTATION = documentation/OptimizingTheMaths.pdf

DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath}))) extern

# PYTHONMAJOR = $(eval python --version | sed -e "s|Python \([0-9]*\)\.[0-9]*\.[0-9]*|\1|g")
# PYTHONMINOR = $(eval python --version | sed -e "s|Python [0-9]*\.\([0-9]*\)\.[0-9]*|\1|g")

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
# 	@echo "  - make deps           - Build a local copy of the external dependencies"

cpp: ${EXECUTABLES} ${PYTHON_LIBRARY_FILE}
doxygen: ${DOXYGEN} 
docs: ${DOCUMENTATION}

# deps: extern/gsl-2.7.1/.libs/libgsl.so extern/boost_1_81_0/stage/lib/libboost_system.so

# Adding the includes and libs for the locally built deps first
FLAGS = -L${CONDA_PREFIX}/lib -Iinclude -I${CONDA_PREFIX}/include -I${CONDA_PREFIX}/include/boost -I${CONDA_PREFIX}/include/${LIBPYTHON}  \
        -lgsl -lgslcblas -l${LIBBOOSTPYTHON} -lboost_program_options -l${LIBPYTHON} -lm -lpthread  \
        -g -std=c++14 -march=native -O3 -MMD -MP \
        \
        -Wno-deprecated-declarations
# Hide the annoying boost auto_ptr=>unique_ptr warning     

ifeq (verbose, $(filter verbose,$(MAKECMDGOALS)))

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	test -n "${CONDA_PREFIX}" || (echo "CONDA_PREFIX not set: Please activate conda environment" ; exit 1)
	g++ $< -o $@ -c ${FLAGS} -fPIC

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	g++ $^ -o $@ ${FLAGS}

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

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	@echo "Building Executable   | g++ ... -o $@"
	@g++ $^ -o $@ ${FLAGS}

endif


${PYTHON_LIBRARY_FILE}: ${LIBRARY_OBJECT_FILES}
	g++ $^ -o $@ -shared ${FLAGS}


${DIRECTORIES}:
	@echo "Making directory      | mkdir -p $@"
	@mkdir -p $@

${DOXYGEN}: ${HEADERS} ${LIBRARY_SOURCES} ${EXECUTABLE_SOURCES}
	@echo "Generating Doxygen Documentation: doxygen Doxyfile ---> $@"
	@doxygen utilities/Doxyfile
	@make -C .doxygen/latex
	@cp .doxygen/latex/refman.pdf $@

#	@python TexLinker.py "https://github.com/Cefhalic/BayesianClusters/blob/" `git rev-parse --abbrev-ref HEAD`


${DOCUMENTATION}:
	@echo "Generating Maths Documentation: pdflatex ... documentation/OptimizingTheMaths ---> documentation/OptimizingTheMaths.pdf"
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths



# extern/gsl-2.7.1.tar.gz: extern
# 	wget -nc https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz -P extern

# extern/gsl-2.7.1/Makefile : extern/gsl-2.7.1.tar.gz
# 	cd extern; \
# 	gtar xzf gsl-2.7.1.tar.gz --skip-old-files

# extern/gsl-2.7.1/.libs/libgsl.so: extern/gsl-2.7.1/Makefile
# 	cd extern/gsl-2.7.1; \
# 	./configure; \
# 	make -j8


# extern/boost_1_81_0.tar.gz: extern
# 	wget -nc https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz -P extern

# extern/boost_1_81_0/bootstrap.sh: extern/boost_1_81_0.tar.gz
# 	cd extern; \
# 	gtar xzf boost_1_81_0.tar.gz --skip-old-files

# extern/boost_1_81_0/stage/lib/libboost_system.so: extern/boost_1_81_0/bootstrap.sh
# 	cd extern/boost_1_81_0; \
# 	./bootstrap.sh; \
# 	./b2 --with-python --with-system --with-program_options 
