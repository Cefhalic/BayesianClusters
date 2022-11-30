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

DIRECTORIES = $(sort $(foreach filePath,${LIBRARY_OBJECT_FILES} ${EXECUTABLE_OBJECT_FILES}, $(dir ${filePath})))

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

FLAGS = -I/usr/include -I/usr/include/python3.9 \
        -L/usr/lib -L/usr/lib64 \
        -lgsl -lgslcblas -lboost_python39 -lpython3.9 \
        -g -std=c++11 -march=native -O3 -lm -lgsl -lpthread\
        -MMD -MP -lboost_program_options

ifeq (verbose, $(filter verbose,$(MAKECMDGOALS)))

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	g++ $< -o $@ -c ${FLAGS} -Iinclude -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	g++ $< -o $@ -c ${FLAGS} -Iinclude -fPIC

-include $(LIBRARY_OBJECT_FILES:.o=.d)
-include $(EXECUTABLE_OBJECT_FILES:.o=.d)

${EXECUTABLES}: %.exe: obj/bin/%.o ${LIBRARY_OBJECT_FILES}
	g++ $^ -o $@ ${FLAGS}

else

.SECONDEXPANSION:
obj/bin/%.o : src/%.cxx | $$(dir obj/bin/%.o)
	@echo "Building Object Files | g++ -c ... $< -o $@"
	@g++ $< -o $@ -c ${FLAGS} -Iinclude -fPIC

.SECONDEXPANSION:
obj/lib/%.o : src/%.cpp | $$(dir obj/lib/%.o)
	@echo "Building Object Files | g++ -c ... $< -o $@"
	@g++ $< -o $@ -c ${FLAGS} -Iinclude -fPIC

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
