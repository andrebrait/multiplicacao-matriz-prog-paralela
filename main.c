#include <stdio.h>
#include <stdlib.h>
#include "dgemm-blocked.h"

/*
Análise:

sudo perf record -e L1-dcache-load-misses main/main
sudo perf report -v
*/

#define LEN 100

void print_matriz(int n, double** matriz);
void fill(int n, double** matriz);
double** square_malloc(int n);

double **A, **B, **C;

int main(int argc, char** argv){
    
    A = square_malloc(LEN);
    B = square_malloc(LEN);
    C = square_malloc(LEN);

    fill(LEN, A);
    fill(LEN, B);
    square_dgemm(LEN, A, B, C);
    print_matriz(LEN, C);
    
    return 0;
}

double** square_malloc(int n){
    double **mat = malloc(n * sizeof(double*));
    int i;
    for(i = 0; i < n; i++){
        mat[i] = malloc(n * sizeof(double));
    }
    return mat;
}

void fill(int n, double **matriz){
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            matriz[i][j] = (double) (i + j);
        }
    }
}

void print_matriz(int n, double **matriz){
    int i, j;
    printf("{\n");
    for(i = 0; i < n; i++){
        printf("{");
        for(j = 0; j < n; j++){
            printf("%f", matriz[i][j]);
            if(j < LEN - 1){
                printf(", ");
            } else {
                printf("}");
            }
        }
        if(i < LEN - 1){
            printf(",\n");
        } else {
            printf("\n}");
        }
    }
    printf("\n");
}
