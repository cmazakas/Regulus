MKFILE    = Makefile

GCC       = g++ -Ofast -march=native -Wall -Wextra -std=c++11 -pedantic -lpthread -lgmp -lmpfr `pkg-config --cflags --libs eigen3`

CSOURCE   = main.cpp point.cpp peano.cpp tetra.cpp neighbour.cpp delaunay.cpp
CHEADER   = structures.hpp mpreal.h
OBJECTS   = ${CSOURCE:.cpp=.o}
EXECBIN   = regulus

all : ${EXECBIN}

${EXECBIN} : ${OBJECTS}
	${GCC} -o $@ ${OBJECTS} -lm

 
%.o : %.cpp

	${GCC} -c $<


clean :
	- rm ${OBJECTS} ${EXECBIN}


again :
	${MAKE} clean all


