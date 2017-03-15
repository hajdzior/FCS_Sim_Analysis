PROG = simFCS
CC = g++
CFLAGS_WARNINGS  = -Wall -Wextra -ansi -pedantic
CFLAGS_DEBUG  = -g
CFLAGS_RELEASE  = -O3 -march=native
LIBS   = -lm -lgsl -lgslcblas

ifeq ($(MAKECMDGOALS),debug)
CFLAGS = $(CFLAGS_DEBUG) $(CFLAGS_WARNINGS)
else
CFLAGS = $(CFLAGS_WARNINGS) $(CFLAGS_RELEASE)
endif

.PHONY: debug release clean run

release: $(PROG)

debug: $(PROG)

$(PROG): simFCS.o 
	${CC} -o $(PROG) simFCS.o ${LIBS}

simFCS.o: simFCS.cpp config.hpp FCS.hpp
	${CC} $(CFLAGS) -c simFCS.cpp -o simFCS.o 

clean:
	rm *.o $(PROG)

run:	$(PROG)
	./$(PROG)
	
