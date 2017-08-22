#include "dgemm-blocked.h"

void clean(int n, double **A){
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            A[i][j] = 0.0;
        }
    }
}

#ifndef slow
void square_dgemm(int n, double **A, double **B, double **C){
    clean(n, C);
    int i, j, k; 
    double r;
    for(i = 0; i < n; i++){
        for (k = 0; k < n; k++){
            r = A[i][k];
            for (j = 0; j < n; j++){
                C[i][j] += r * B[k][j];
            }
        }
    }
}
#else
void square_dgemm(int n, double **A, double **B, double **C){
    clean(n, C);
    int i, j, k; 
    double r;
    for(j = 0; j < n; j++){
        for (k = 0; k < n; k++){
            r = B[k][j];
            for (i = 0; i < n; i++){
                C[i][j] += r * A[i][k];
            }
        }
    }
}
#endif