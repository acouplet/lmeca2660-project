CC=gcc
CFLAGS=-O3 -W -Wall -lm -g
LDFLAGS=
EXEC=main

all: $(EXEC)

main: cfd.o Problem.o main.o
	$(CC) -o main cfd.o Problem.o main.o $(LDFLAGS) $(CFLAGS)

cfd.o: cfd.c
	$(CC) -o cfd.o -c cfd.c $(CFLAGS)

Problem.o: Problem.c cfd.h
	$(CC) -o Problem.o -c Problem.c $(CFLAGS)

main.o: main.c Problem.h cfd.h
	$(CC) -o main.o -c main.c $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
