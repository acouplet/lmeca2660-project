CC=gcc
CFLAGS=-O0 -W -Wall -lm -g
LDFLAGS=
EXEC=main

all: $(EXEC)

main: cfd.o main.o
	$(CC) -o main cfd.o main.o $(LDFLAGS) $(CFLAGS)

cfd.o: cfd.c
	$(CC) -o cfd.o -c cfd.c $(CFLAGS)

main.o: main.c cfd.h
	$(CC) -o main.o -c main.c $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
