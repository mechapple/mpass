CC = g++
CPP_FLAGS = -std=c++0x

all: main debug

main: bernstein.h  bezier.h  converged.h  cubature.h  firemin.h  lj.h  point.h  vwrapper.h
	$(CC) $(CPP_FLAGS) -o main main.cxx hcubature.c pcubature.c

debug: bernstein.h  bezier.h  converged.h  cubature.h  firemin.h  lj.h  point.h  vwrapper.h
	$(CC) $(CPP_FLAGS) -o debug debug.cxx hcubature.c pcubature.c
