#to compile the executable for dedup, the redundant contig remover

CC = g++
CFLAGS = -g -Wall -std=c++0x 

default: dedup
dedup: dlib.o dedup.o
	$(CC) $(CFLAGS) -o dedup dlib.o dedup.o

dlib.o: dlib.cpp dedup.h
	$(CC) $(CFLAGS) -c dlib.cpp

dedup.o: dedup.cpp dedup.h seqIO.h
	$(CC) $(CFLAGS) -c dedup.cpp

clean:
	$(RM) *.o 
