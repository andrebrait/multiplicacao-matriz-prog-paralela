#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <immintrin.h>

const char *dgemm_desc = "Vectorial Matmult";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 64
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))

#include <emmintrin.h>

#ifdef USE_RESTRICT
#else
#define restrict
#endif

/*
 * Calculates the transpose of a matrix M.
*/
double *transpose(int lda, double *restrict M) {
	double *result = (double *) malloc(sizeof(double) * lda * lda);
	for(int i = 0; i < lda; i++)
		for(int j = 0; j < lda; j++)
			result[j * lda + i] = M[i * lda + j];
	return result;
}

/* 
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. 
 * 1st. loop: for each block-row of A;
 * 2nd. loop: for each block-column of B;
 * 3rd. loop: accumulate block dgemms into block of C.
 *   Then, we correct the block dimensions if block "goes off edge of" the matrix
 *   and perform an individual block dgemm.
 */
void square_dgemm(int lda, double *restrict A, double *restrict B, double *restrict C) {
	B = transpose(lda, B);
	for(int i = 0; i < lda; i++) {
		
	}
}
