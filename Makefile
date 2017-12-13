# on Edision we will benchmark you against the default vendor-tuned BLAS. The compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. To change compilers, you need to type something like: "module swap PrgEnv-pgi PrgEnv-gnu" See the NERSC documentation for available compilers.

CC = cc 
OPT = -O0
CFLAGS = -Wall -std=gnu99 -fstrict-aliasing $(OPT)
LDFLAGS = -Wall
# librt is needed for clock_gettime
LDLIBS = -lrt -lblas -lm -msse2 -msse3 -msse4 -mavx

targets = benchmark-naive benchmark-blocked benchmark-blas benchmark-strassen benchmark-vect
objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o  dgemm-strassen.o dgemm-vect.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-strassen : benchmark.o dgemm-strassen.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-vect: benchmark.o dgemm-vect.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
