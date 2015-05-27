CFLAGS=-ansi -std=c99 -W -Wall -O2 $(DEFINEA) $(DEFINEB)
CCFLAGS=-ansi -W -Wall -O2 $(DEFINEA) $(DEFINEB)
CCC=g++
CC=gcc
EXE=

OBJ=mathtext.o gvars.o gmaths.o erime.o rime.o
LIB=librime.a
TEST=test

lib: $(OBJ)
	ar rcvs $(LIB) $(OBJ)

.cc.o:
	$(CCC) $(CCFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

test: $(LIB)
	$(CCC) $(CCFLAGS) $(TEST).cc -L. -lrime -L/usr/lib -lgsl -lgslcblas -lm -o $(TEST)$(EXE)

clean:
	rm $(OBJ)

mathtext.o: ./mathtext.h
gvars.o: ./gvars.h
gmaths.o: ./gvars.h ./gmaths.h
erime.o: ./gvars.h ./gmaths.h ./mathtext.h
rime.o: ./erime.h ./gvars.h ./gmaths.h ./mathtext.h
librime.a: ./rime.h ./erime.h ./gvars.h ./gmaths.h ./mathtext.h
