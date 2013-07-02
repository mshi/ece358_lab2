CXX           = g++					# compiler
CXXFLAGS      = -g -Wno-unused-label -MMD
MAKEFILE_NAME = ${firstword ${MAKEFILE_LIST}}	# makefile name
OBJECTS       = lab2.o random.o Debug.o simulation.o
DEPENDS 			= ${OBJECTS:.o=.d}			# substitute ".o" with ".d"
EXEC          = lab2

.PHONY: clean

all: ${EXEC}

${EXEC}: ${OBJECTS}
	${CXX} $^ -o $@

${OBJECTS} : ${MAKEFILE_NAME}			# OPTIONAL : changes to this file => recompile

-include ${DEPENDS}

clean:
	rm -f *.d *.o ${EXEC}
