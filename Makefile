CC = mpicc
CFLAGS = -g -lm -Wall -std=c99
OBJECTS = multiple.o
all : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o multiple


multiple: .c
	$(CC) $(CFLAGS) -c multiple.c

.PHONY : clean
clean :
	-rm edit $(OBJECTS)
