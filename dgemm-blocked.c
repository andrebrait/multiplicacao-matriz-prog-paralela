#include "dgemm-blocked.h"

void clean(int, double **);
void square_dgemm_ikj(int, double **, double **, double **);
void square_dgemm_jki(int, double **, double **, double **);

void square_dgemm(int n, double **A, double **B, double **C) {
#ifdef IKJ
    square_dgemm_ikj(n, A, B, C);
#elif JKI
    square_dgemm_jki(n, A, B, C);
#endif
}

void clean(int n, double **A) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = 0.0;
        }
    }
}

void square_dgemm_ikj(int n, double **A, double **B, double **C) {
    clean(n, C);
    int i, j, k;
    double r;
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            r = A[i][k];
            for (j = 0; j < n; j++) {
                C[i][j] += r * B[k][j];
            }
        }
    }
}

void square_dgemm_jki(int n, double **A, double **B, double **C) {
    clean(n, C);
    int i, j, k;
    double r;
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
            r = B[k][j];
            for (i = 0; i < n; i++) {
                C[i][j] += r * A[i][k];
            }
        }
    }
}
