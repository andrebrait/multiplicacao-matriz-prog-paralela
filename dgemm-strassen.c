#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define SubM_11_cm(i, j, _n, n) (j * n + i)
#define SubM_12_cm(i, j, _n, n) ((j + _n) * n + i)
#define SubM_21_cm(i, j, _n, n) (j * n + (i + _n))
#define SubM_22_cm(i, j, _n, n) ((j + _n) * n + (i + _n))

const char *dgemm_desc = "Strassen divide and conquer dgemm.";

/**
* Function signatures.
*/
static void dgemm_strassen(double *restrict A, double *restrict B, double *restrict C, int n);

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
	double *restrict M1, *restrict M2, *restrict M3, *restrict M4, *restrict M5, *restrict M6, *restrict M7;
	
	double *restrict A_11 = createMatrixColumnMajor(_n);
	double *restrict A_22 = createMatrixColumnMajor(_n);
	double *restrict B_11 = createMatrixColumnMajor(_n);
	double *restrict B_22 = createMatrixColumnMajor(_n);
	double *restrict T0 = createMatrixColumnMajor(_n);
	double *restrict T1 = createMatrixColumnMajor(_n);
	double *restrict T2 = createMatrixColumnMajor(_n);
	double *restrict T3 = createMatrixColumnMajor(_n);
	double *restrict T4 = createMatrixColumnMajor(_n);
	double *restrict T5 = createMatrixColumnMajor(_n);
	double *restrict T6 = createMatrixColumnMajor(_n);
	double *restrict T7 = createMatrixColumnMajor(_n);
	double *restrict T8 = createMatrixColumnMajor(_n);
	double *restrict T9 = createMatrixColumnMajor(_n);

	
	M1 = createMatrixColumnMajor(_n);
	M2 = createMatrixColumnMajor(_n);
	M3 = createMatrixColumnMajor(_n);
	M4 = createMatrixColumnMajor(_n);
	M5 = createMatrixColumnMajor(_n);
	M6 = createMatrixColumnMajor(_n);
	M7 = createMatrixColumnMajor(_n);
	
	register int i, j;
	for(j = 0; j < _n; j++) {
		for(i = 0; i < _n; i++) {
			A_11[j * _n + i] = A[SubM_11_cm(i, j, _n, n)];
			A_22[j * _n + i] = A[SubM_22_cm(i, j, _n, n)];
			B_11[j * _n + i] = B[SubM_11_cm(i, j, _n, n)];
			B_22[j * _n + i] = B[SubM_22_cm(i, j, _n, n)];
			T1[j * _n + i] = A[SubM_11_cm(i, j, _n, n)] + A[SubM_22_cm(i, j, _n, n)];
			T2[j * _n + i] = B[SubM_11_cm(i, j, _n, n)] + B[SubM_22_cm(i, j, _n, n)];
			T3[j * _n + i] = A[SubM_21_cm(i, j, _n, n)] + A[SubM_22_cm(i, j, _n, n)];
			T4[j * _n + i] = B[SubM_12_cm(i, j, _n, n)] - B[SubM_22_cm(i, j, _n, n)];
			T5[j * _n + i] = B[SubM_21_cm(i, j, _n, n)] - B[SubM_11_cm(i, j, _n, n)];
			T6[j * _n + i] = A[SubM_11_cm(i, j, _n, n)] + A[SubM_12_cm(i, j, _n, n)];
			T7[j * _n + i] = A[SubM_21_cm(i, j, _n, n)] - A[SubM_11_cm(i, j, _n, n)];
			T8[j * _n + i] = B[SubM_11_cm(i, j, _n, n)] + B[SubM_12_cm(i, j, _n, n)];
			T9[j * _n + i] = A[SubM_12_cm(i, j, _n, n)] - A[SubM_22_cm(i, j, _n, n)];
			T0[j * _n + i] = B[SubM_21_cm(i, j, _n, n)] + B[SubM_22_cm(i, j, _n, n)];
		}
	}
	
	dgemm_strassen(T1, T2, M1, _n);
	dgemm_strassen(T3, B_11, M2, _n);	
	dgemm_strassen(A_11, T4, M3, _n);
	dgemm_strassen(A_22, T5, M4, _n);
	dgemm_strassen(T6, B_22, M5, _n);
	dgemm_strassen(T7, T8 , M6, _n);
	dgemm_strassen(T9, T0, M7, _n);
	
	for(j = 0; j < _n; j++) {
		for(i = 0; i < _n; i++) {
			C[SubM_11_cm(i, j, _n, n)] = M1[j * _n + i] + M4[j * _n + i] - M5[j * _n + i] + M7[j * _n + i];
			C[SubM_12_cm(i, j, _n, n)] = M3[j * _n + i] + M5[j * _n + i];
			C[SubM_21_cm(i, j, _n, n)] = M2[j * _n + i] + M4[j * _n + i];
			C[SubM_22_cm(i, j, _n, n)] = M1[j * _n + i] - M2[j * _n + i] + M3[j * _n + i] + M6[j * _n + i];
		}
	}

	A_11 = freeMatrixColumnMajor(A_11, _n);
	A_22 = freeMatrixColumnMajor(A_22, _n);
	B_11 = freeMatrixColumnMajor(B_11, _n);
	B_22 = freeMatrixColumnMajor(B_22, _n);
	M1 = freeMatrixColumnMajor(M1, _n);
	M2 = freeMatrixColumnMajor(M2, _n);
	M3 = freeMatrixColumnMajor(M3, _n);
	M4 = freeMatrixColumnMajor(M4, _n);
	M5 = freeMatrixColumnMajor(M5, _n);
	M6 = freeMatrixColumnMajor(M6, _n);
	M7 = freeMatrixColumnMajor(M7, _n);
	T0 = freeMatrixColumnMajor(T0, _n);
	T1 = freeMatrixColumnMajor(T1, _n);
	T2 = freeMatrixColumnMajor(T2, _n);
	T3 = freeMatrixColumnMajor(T3, _n);
	T4 = freeMatrixColumnMajor(T4, _n);
	T5 = freeMatrixColumnMajor(T5, _n);
	T6 = freeMatrixColumnMajor(T6, _n);
	T7 = freeMatrixColumnMajor(T7, _n);
	T8 = freeMatrixColumnMajor(T8, _n);
	T9 = freeMatrixColumnMajor(T9, _n);
}

/**
* Creates a real matrix size * size in column-major format.
*/
static double *restrict createMatrixColumnMajor(int size) {
	double *restrict M = (double *) malloc(sizeof(double) * size * size);
	if(M == NULL) {
		exit(1);
	}
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
	for(j = 0; j <  correctSize; j++) {
		for(i = 0; i < correctSize; i++) {
			if(j >= n || i >= n) {
				Asized2[j * correctSize + i] = 0.0;
				Bsized2[j * correctSize + i] = 0.0;			
			}
			else {
				Asized2[j * correctSize + i] = A[j * n + i];
				Bsized2[j * correctSize + i] = B[j * n + i];
			}
		}
	}
	double *restrict Csized2 = createMatrixColumnMajor(correctSize);
	dgemm_strassen(Asized2, Bsized2, Csized2, correctSize);
	for(j = 0; j < correctSize; j++) {
		for(i = 0; i < correctSize; i++) {
			if(i < n && j < n) {
				C[j * n + i] = Csized2[j * correctSize + i];				
			}
		}
	}
	freeMatrixColumnMajor(Asized2, correctSize);
	freeMatrixColumnMajor(Bsized2, correctSize);
	freeMatrixColumnMajor(Csized2, correctSize);
}
