#include <math.h>
#include <stdlib.h>
#include <string.h>

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

const char *dgemm_desc = "Strassen divide and conquer dgemm.";

/**
* Function signatures.
*/
static void dgemm_strassen(double **A, double **B, double **C, int n);

static void plus(double **A, double **B, double **A_plus_B, int size);

static void minus(double **A, double **B, double **A_minus_B, int size);

static double **createMatrix(int size);

static double **freeMatrix(double **M, int size);

/**
* Returns the first value 2^k >= n.
*/
static int next_power_of_two(int n) {
    return (int)(pow(2, ceil(log((double)n) / log(2.0))));
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
void dgemm_strassen(double **A, double **B, double **C, int n) {
	if(n == 1) {
		C[0][0] = A[0][0] * B[0][0];
		return;
	}
	int _n = n/2;
	double **A_11, **A_12, **A_21, **A_22;
	double **B_11, **B_12, **B_21, **B_22;
	double **C_11, **C_12, **C_21, **C_22;
	double **M1, **M2, **M3, **M4, **M5, **M6, **M7;
	
	A_11 = createMatrix(_n);
	A_12 = createMatrix(_n);
	A_21 = createMatrix(_n);
	A_22 = createMatrix(_n);
	
	B_11 = createMatrix(_n);
	B_12 = createMatrix(_n);
	B_21 = createMatrix(_n);
	B_22 = createMatrix(_n);
	
	C_11 = createMatrix(_n);
	C_12 = createMatrix(_n);
	C_21 = createMatrix(_n);
	C_22 = createMatrix(_n);
	
	M1 = createMatrix(_n);
	M2 = createMatrix(_n);
	M3 = createMatrix(_n);
	M4 = createMatrix(_n);
	M5 = createMatrix(_n);
	M6 = createMatrix(_n);
	M7 = createMatrix(_n);
	
	double **APartialResult = createMatrix(_n);
	double **BPartialResult = createMatrix(_n);
	
	register int i, j;
	
	for(i = 0; i < _n; i++) {
		for(j = 0; j < _n; j++) {
			A_11[i][j] = A[i][j];
			A_12[i][j] = A[i][j + _n];
			A_21[i][j] = A[i + _n][j];
			A_22[i][j] = A[i + _n][j + _n];
			
			B_11[i][j] = B[i][j];
			B_12[i][j] = B[i][j + _n];
			B_21[i][j] = B[i + _n][j];
			B_22[i][j] = B[i + _n][j + _n];
		}
	}
	
	plus(A_11, A_22, APartialResult, _n);
	plus(B_11, B_22, BPartialResult, _n);
	dgemm_strassen(APartialResult, BPartialResult, M1, _n);
	
	plus(A_21, A_22, APartialResult, _n);
	dgemm_strassen(APartialResult, B_11, M2, _n);
	
	minus(B_12, B_22, BPartialResult, _n);
	dgemm_strassen(A_22, BPartialResult, M3, _n);
	
	minus(B_21, B_11, BPartialResult, _n);
	dgemm_strassen(A_22, BPartialResult, M4, _n);
	
	plus(A_11, A_12, APartialResult, _n);
	dgemm_strassen(APartialResult, B_22, M5, _n);
	
	minus(A_21, A_11, APartialResult, _n);
	plus(B_11, B_12, BPartialResult, _n);
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
			C[i][j] = C_11[i][j];
			C[i][j + _n] = C_12[i][j];
			C[i + _n][j] = C_21[i][j];
			C[i + _n][j + _n] = C_22[i][j];
		}
	}
	
	A_11 = freeMatrix(A_11, _n);
	A_12 = freeMatrix(A_12, _n);
	A_21 = freeMatrix(A_21, _n);
	A_22 = freeMatrix(A_22, _n);
	
	B_11 = freeMatrix(B_11, _n);
	B_12 = freeMatrix(B_12, _n);
	B_21 = freeMatrix(B_21, _n);
	B_22 = freeMatrix(B_22, _n);
	
	C_11 = freeMatrix(C_11, _n);
	C_12 = freeMatrix(C_12, _n);
	C_21 = freeMatrix(C_21, _n);
	C_22 = freeMatrix(C_22, _n);
	
	M1 = freeMatrix(M1, _n);
	M2 = freeMatrix(M2, _n);
	M3 = freeMatrix(M3, _n);
	M4 = freeMatrix(M4, _n);
	M5 = freeMatrix(M5, _n);
	M6 = freeMatrix(M6, _n);
	M7 = freeMatrix(M7, _n);
	
	APartialResult = freeMatrix(APartialResult, _n);
	BPartialResult = freeMatrix(BPartialResult, _n);
}

/**
* Given two real matrices A and B, returns A + B where A, B and 
* C are size x size.
*/
void plus(double **A, double **B, double **A_plus_B, int size) {	
	register int i, j;
	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			A_plus_B[i][j] = A[i][j] + B[i][j];
		}
	}
}

/**
* Given two real matrices A and B, returns A - B where A, B and
* C are size x size.
*/
void minus(double **A, double **B, double **A_minus_B, int size) {
	register int i, j;
	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			A_minus_B[i][j] = A[i][j] - B[i][j];	
		}
	}
}

/**
* Creates a real matrix size x size. If there's not enough memory to be
* allocated, the program should be terminated. Otherwise, returns a
* pointer to the created matrix's memory address.
*/
double **createMatrix(int size) {
	double **M = (double **) malloc(sizeof(double *) * size);
	if(M == NULL) {
		exit(1);
	}
	register int i;
	for(i = 0; i < size; i++) {
		M[i] = (double *) malloc(sizeof(double) * size);
		if(M[i] == NULL) {
			exit(1);
		}
		memset(M[i], 0.0, sizeof(double) * size);
	}
	return M;
}

/**
* Frees the allocated space for a real matrix addressed by the
* pointer M, whose size must be size x size. Returns NULL in
* case of success.
*/
double **freeMatrix(double **M, int size) {
	register int i;
	if(M == NULL)
		return NULL;
	for(i = 0; i < size; i++) {
		if(M[i]) {
			free(M[i]);
			M[i] = NULL;
		}
	}
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
	double **Asized2;
	double **Bsized2;
	if(correctSize != n) {
		Asized2 = createMatrix(correctSize);
		Bsized2 = createMatrix(correctSize);
		register int i, j;
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				Asized2[i][j] = A[i][j];
				Bsized2[i][j] = B[i][j];
			}
		}
	}
	else {
		Asized2 = A;
		Bsized2 = B;
	}
	double **Csized2 = createMatrix(correctSize);
	dgemm_strassen(Asized2, Bsized2, Csized2, correctSize);
	if(correctSize != n) {
		register int i, j;
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				C[i][j] = Csized2[i][j];
			}
		}
	}
}
