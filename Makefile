CC            = gcc
CFLAGS       = -ansi -Wall -Wextra -Wno-unused-parameter
LIBS         = -lm

# ===================================================================== #
# ! - all : Compile all the executables
# ! - help : Display this help
# ! - clean : remove all the *.o and executables
# ===================================================================== #

all: lma #for Lavenberg Marquardt algorithm

help:
	@grep -E "^# !" Makefile | sed -e 's/# ! -/-/g'

clean:
	rm -f *.o lma

lma: main.o lm_func.o gaussj.o fit_funcs.o  nrutil.o
	$(CC) $(CFLAGS)  -o lma  main.o lm_func.o gaussj.o fit_funcs.o  nrutil.o    $(LIBS)


lm_func.o: lm_func.c nrutil.h gaussj.h 
	$(CC) $(CFLAGS) -c lm_func.c 

main.o: main.c lm_func.h fit_funcs.h  gaussj.h 
	$(CC) $(CFLAGS) -c main.c 

fit_funcs.o: fit_funcs.c 
	$(CC) $(CFLAGS) -c fit_funcs.c

gaussj.o: gaussj.c  nrutil.h
	$(CC) $(CFLAGS) -c gaussj.c

nrutil.o: nrutil.c
	$(CC) $(CFLAGS) -c nrutil.c  





