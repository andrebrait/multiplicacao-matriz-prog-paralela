#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

const char *dgemm_desc = "Strassen divide and conquer dgemm.";

/**
* Function signatures.
*/
static void dgemm_strassen(double *restrict A, double *restrict B, double *restrict C, int n);

static void plus(double *restrict A, double *restrict B, double *restrict A_plus_B, int size);

static void minus(double *restrict A, double *restrict B, double *restrict A_minus_B, int size);

static double *restrict createMatrixColumnMajor(int size);

static double *restrict freeMatrixColumnMajor(double *restrict M, int size);

/**
* Returns the first value 2^k >= n.
*/
static int next_power_of_two(int n) {
    return (int) (pow(2, ceil(log2( (double)n))));
}

/**
* Implementation of Strassen Algorithm for Matrices Multiplication. 
* This algorithm reduces the number of multiplications needed to
* find the result C = AB, where A, B and C are square matrices of
* dimension n. Consider that necessarily exists k such that 2^k = n.
* Then, the algorithm says that
*
*	dgemm_strassen(A, B, C, n) = (C_(0,0) = A_(0,0) * B_(0,0)),
*
* if n = 1. Otherwise, we must rely on that
*
*	C_(1,1) = M1 + M4 - M5 + M7
* 	C_(1,2) = M3 + M5
* 	C_(2,1) = M2 + M4
*	C_(2,2) = M1 - M2 + M3 + M6,
*
* where
*
*	M1 = (A_(1,1) + A_(2,2))*(B_(1,1) + B_(2,2))
*	M2 = (A_(2,1) + A_(2,2))*B_(1,1)
*	M3 = A_(1,1)*(B_(1,2)-B_(2,2))
*	M4 = A_(2,2)*(B_(2,1)-B_(1,1))
*	M5 = (A_(1,1) + A_(1,2))*B_(2,2)
*	M6 = (A_(2,1) - A_(1,1))*(B_(1,1) + B_(1,2))
*	M7 = (A_(2,1) - A_(1,1))*(B_(2,1) + B_(2,2)),
*
* such that every multiplication found in M1, M2, M3, ... are
* resolved by the general call of dgemm_strassen. Note that
* the descrpiption D_(1, .) references a submatrix D(1:n/2,.) and
* that D_(2, .) references a submatrix D(n/2:n,.), where the '.' (dot)
* is equals to "whatever". Therefore, all multiplications in
* dgemm_strassen are the result of a divide and conquer based on
* multiple submatrices of A and B, then partial results of M1, M2, M3, ...
* are computed in each recursive call, so that C is a combination 
* of these results. In other words, the general call for dgemm_srassen
* is
*
*	dgemm_strassen(A, B, C, n) = 
*
*		C_(1,1) = M1 + M4 - M5 + M7
* 		C_(1,2) = M3 + M5
* 		C_(2,1) = M2 + M4
*		C_(2,2) = M1 - M2 + M3 + M6,
*
* 	where each M1, M2, ... are previously resolved in a recursive call of
*
*		dgemm_strassen(A_(1,1) + A_(2,2), B_(1,1) + B_(2,2), M1, n/2)...
*
* ...and all corresponding expression to M2, M3, ... in order to compute
* each submatrix of C.
*/
static void dgemm_strassen(double *restrict A, double *restrict B, double *restrict C, int n) {
	if(n == 1) {
		C[0] = A[0] * B[0];
		return;
	}
	int _n = n/2;
	double *restrict A_11, *restrict A_12, *restrict A_21, *restrict A_22;
	double *restrict B_11, *restrict B_12, *restrict B_21, *restrict B_22;
	double *restrict C_11, *restrict C_12, *restrict C_21, *restrict C_22;
	double *restrict M1, *restrict M2, *restrict M3, *restrict M4, *restrict M5, *restrict M6, *restrict M7;
	
	A_11 = createMatrixColumnMajor(_n);
	A_12 = createMatrixColumnMajor(_n);
	A_21 = createMatrixColumnMajor(_n);
	A_22 = createMatrixColumnMajor(_n);
	
	B_11 = createMatrixColumnMajor(_n);
	B_12 = createMatrixColumnMajor(_n);
	B_21 = createMatrixColumnMajor(_n);
	B_22 = createMatrixColumnMajor(_n);
	
	C_11 = createMatrixColumnMajor(_n);
	C_12 = createMatrixColumnMajor(_n);
	C_21 = createMatrixColumnMajor(_n);
	C_22 = createMatrixColumnMajor(_n);
	
	M1 = createMatrixColumnMajor(_n);
	M2 = createMatrixColumnMajor(_n);
	M3 = createMatrixColumnMajor(_n);
	M4 = createMatrixColumnMajor(_n);
	M5 = createMatrixColumnMajor(_n);
	M6 = createMatrixColumnMajor(_n);
	M7 = createMatrixColumnMajor(_n);
	
	double *restrict APartialResult = createMatrixColumnMajor(_n);
	double *restrict BPartialResult = createMatrixColumnMajor(_n);
	
	register int i, j;
	
	for(i = 0; i < _n; i++) {
		for(j = 0; j < _n; j++) {
			A_11[j * _n + i] = A[j * n + i];
			A_12[j * _n + i] = A[(j + _n) * n + i];
			A_21[j * _n + i] = A[j * n + (i + _n)];
			A_22[j * _n + i] = A[(j + _n) * n + (i + _n)];
			
			B_11[j * _n + i] = B[j * n + i];
			B_12[j * _n + i] = B[(j + _n) * n + i];
			B_21[j * _n + i] = B[j * n + (i + _n)];
			B_22[j * _n + i] = B[(j + _n) * n + (i + _n)];
		}
	}
	
	plus(A_11, A_22, APartialResult, _n);
	plus(B_11, B_22, BPartialResult, _n);
	dgemm_strassen(APartialResult, BPartialResult, M1, _n);
	
	plus(A_21, A_22, APartialResult, _n);
	dgemm_strassen(APartialResult, B_11, M2, _n);
	
	minus(B_12, B_22, BPartialResult, _n);
	dgemm_strassen(A_11, BPartialResult, M3, _n);
	
	minus(B_21, B_11, BPartialResult, _n);
	dgemm_strassen(A_22, BPartialResult, M4, _n);
	
	plus(A_11, A_12, APartialResult, _n);
	dgemm_strassen(APartialResult, B_22, M5, _n);
	
	plus(B_11, B_12, BPartialResult, _n);
	minus(A_21, A_11, APartialResult, _n);
	dgemm_strassen(APartialResult, BPartialResult, M6, _n);
	
	minus(A_12, A_22, APartialResult, _n);
	plus(B_21, B_22, BPartialResult, _n);
	dgemm_strassen(APartialResult, BPartialResult, M7, _n);
	
	plus(M3, M5, C_12, _n);
	plus(M2, M4, C_21, _n);
	
	plus(M1, M4, APartialResult, _n);
	plus(APartialResult, M7, BPartialResult, _n);
	minus(BPartialResult, M5, C_11, _n);
	
	plus(M1, M3, APartialResult, _n);
	plus(APartialResult, M6, BPartialResult, _n);
	minus(BPartialResult, M2, C_22, _n);
	
	for(i = 0; i < _n; i++) {
		for(j = 0; j < _n; j++) {
			C[j * n + i] = C_11[j * _n + i];
			C[(j + _n) * n + i] = C_12[j * _n + i];
			C[j * n + (i + _n)] = C_21[j * _n + i];
			C[(j + _n) * n + (i + _n)] = C_22[j * _n + i];
		}
	}
	
	A_11 = freeMatrixColumnMajor(A_11, _n);
	A_12 = freeMatrixColumnMajor(A_12, _n);
	A_21 = freeMatrixColumnMajor(A_21, _n);
	A_22 = freeMatrixColumnMajor(A_22, _n);
	
	B_11 = freeMatrixColumnMajor(B_11, _n);
	B_12 = freeMatrixColumnMajor(B_12, _n);
	B_21 = freeMatrixColumnMajor(B_21, _n);
	B_22 = freeMatrixColumnMajor(B_22, _n);
	
	C_11 = freeMatrixColumnMajor(C_11, _n);
	C_12 = freeMatrixColumnMajor(C_12, _n);
	C_21 = freeMatrixColumnMajor(C_21, _n);
	C_22 = freeMatrixColumnMajor(C_22, _n);
	
	M1 = freeMatrixColumnMajor(M1, _n);
	M2 = freeMatrixColumnMajor(M2, _n);
	M3 = freeMatrixColumnMajor(M3, _n);
	M4 = freeMatrixColumnMajor(M4, _n);
	M5 = freeMatrixColumnMajor(M5, _n);
	M6 = freeMatrixColumnMajor(M6, _n);
	M7 = freeMatrixColumnMajor(M7, _n);
	
	APartialResult = freeMatrixColumnMajor(APartialResult, _n);
	BPartialResult = freeMatrixColumnMajor(BPartialResult, _n);
}

/**
* Given two real matrices A and B, returns A + B where A, B and 
* C are size x size.
*/
static void plus(double *restrict A, double *restrict B, double *restrict A_plus_B, int size) {	
	register int i, j;
	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			A_plus_B[j * size + i] = A[j * size + i] + B[j * size + i];
		}
	}
}

/**
* Given two real matrices A and B, returns A - B where A, B and
* C are size x size.
*/
static void minus(double *restrict A, double *restrict B, double *restrict A_minus_B, int size) {
	register int i, j;
	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			A_minus_B[j * size + i] = A[j * size + i] - B[j * size + i];	
		}
	}
}

/**
* Creates a real matrix size * size in column-major format.
*/
static double *restrict createMatrixColumnMajor(int size) {
	double *restrict M = (double *) malloc(sizeof(double) * size * size);
	if(M == NULL) {
		exit(1);
	}
	register int i;
	memset(M, 0.0, sizeof(double) * size * size);
	return M;
}

/**
* Frees the allocated space for a real matrix in column-major format.
*/
static double *restrict freeMatrixColumnMajor(double *restrict M, int size) {
	if(M == NULL)
		return NULL;
	free(M);
	M = NULL;
	return NULL;
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, double *restrict A, double *restrict B, double *restrict C) {
	int correctSize = next_power_of_two(n);
	double *restrict Asized2 = createMatrixColumnMajor(correctSize);
	double *restrict Bsized2 = createMatrixColumnMajor(correctSize);
	register int i, j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			Asized2[j * n + i] = A[j * n + i];
			Bsized2[j * n + i] = B[j * n + i];
		}
	}
	double *restrict Csized2 = createMatrixColumnMajor(correctSize);
	dgemm_strassen(Asized2, Bsized2, Csized2, correctSize);
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			C[j * n + i] = Csized2[j * n + i];
		}
	}
	freeMatrixColumnMajor(Asized2, correctSize);
	freeMatrixColumnMajor(Bsized2, correctSize);
	freeMatrixColumnMajor(Csized2, correctSize);
}
