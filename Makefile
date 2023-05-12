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

LIBPYTHON = $(shell python -c "from sys import version_info; print( f'python{version_info[0]}.{version_info[1]}' )" )
LIBBOOSTPYTHON = $(shell python -c "from sys import version_info; print( f'boost_python{version_info[0]}{version_info[1]}' )" )



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
	@echo "  - make deps           - Build a local copy of the external dependencies"
	@echo

cpp: ${EXECUTABLES} ${PYTHON_LIBRARY_FILE}
doxygen: ${DOXYGEN} 
docs: ${DOCUMENTATION}

deps: extern/gsl-2.7.1/.libs/libgsl.so extern/boost_1_81_0/stage/lib/libboost_system.so

# Adding the includes and libs for the locally built deps first
FLAGS = -Iextern/gsl-2.7.1 -Iextern/boost_1_81_0 \
	-Lextern/gsl-2.7.1/cblas/.libs -Lextern/gsl-2.7.1/.libs -Lextern/boost_1_81_0/stage/lib \
	\
	-Iinclude -I/usr/include/${LIBPYTHON} \
	-lgsl -lgslcblas -l${LIBBOOSTPYTHON} -lboost_program_options -l${LIBPYTHON} -lm -lpthread \
        -g -std=c++11 -march=native -O3 -MMD -MP


ifeq (verbose, $(filter verbose,$(MAKECMDGOALS)))

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	g++ $< -o $@ -c ${FLAGS} -fPIC

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	g++ $^ -o $@ ${FLAGS}

else

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	@echo "Building Object Files | g++ -c ... $< -o $@"
	@g++ $< -o $@ -c ${FLAGS} -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
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
	@doxygen Doxyfile
	@make -C .doxygen/latex
	@cp .doxygen/latex/refman.pdf $@

#	@python TexLinker.py "https://github.com/Cefhalic/BayesianClusters/blob/" `git rev-parse --abbrev-ref HEAD`


${DOCUMENTATION}:
	@echo "Generating Maths Documentation: pdflatex ... documentation/OptimizingTheMaths ---> documentation/OptimizingTheMaths.pdf"
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths
	@pdflatex -output-directory=./documentation ./documentation/OptimizingTheMaths



extern/gsl-2.7.1.tar.gz: extern
	wget -nc https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz -P extern

extern/gsl-2.7.1/Makefile : extern/gsl-2.7.1.tar.gz
	cd extern; \
	gtar xzf gsl-2.7.1.tar.gz --skip-old-files

extern/gsl-2.7.1/.libs/libgsl.so: extern/gsl-2.7.1/Makefile
	cd extern/gsl-2.7.1; \
	./configure; \
	make -j8


extern/boost_1_81_0.tar.gz: extern
	wget -nc https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz -P extern

extern/boost_1_81_0/bootstrap.sh: extern/boost_1_81_0.tar.gz
	cd extern; \
	gtar xzf boost_1_81_0.tar.gz --skip-old-files

extern/boost_1_81_0/stage/lib/libboost_system.so: extern/boost_1_81_0/bootstrap.sh
	cd extern/boost_1_81_0; \
	./bootstrap.sh; \
	./b2 --with-python --with-system --with-program_options 
