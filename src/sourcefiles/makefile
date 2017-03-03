CFLAGS=-g -Wall -Wextra -Isrc -rdynamic -DNDEBUG $(OPTFLAGS)

#all:
	#rm -f main
	#make main

all: simplex.o interface.o matrix.o main.o
	gcc -o main simplex.o interface.o matrix.o main.o
	rm simplex.o interface.o matrix.o main.o

simplex.o: matrix.o interface.o
	gcc -c simplex.c

interface.o:
	gcc -c interface.c

matrix.o:
	gcc -c matrix.c
