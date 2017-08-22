#include<stdio.h>

#define LEN 100

/*
Análise:

sudo perf record -e L1-dcache-load-misses ./main
sudo perf report -v
*/

void multiplica_matriz(int matrizA[LEN][LEN], int matrizB[LEN][LEN], int matrizC[LEN][LEN]);
void multiplica_matriz_slow(int matrizA[LEN][LEN], int matrizB[LEN][LEN], int matrizC[LEN][LEN]);
void print_matriz(int matriz[LEN][LEN]);
void fill(int matriz[LEN][LEN]);
void clean(int matriz[LEN][LEN]);

int matrizA[LEN][LEN], matrizB[LEN][LEN], matrizC[LEN][LEN];

int main(int argc, char** argv){
    
    fill(matrizA);
    fill(matrizB);
    
    multiplica_matriz(matrizA, matrizB, matrizC);
    //print_matriz(matrizC);

    return 0;
}

void fill(int matriz[LEN][LEN]){
    int i, j;
    for(i = 0; i < LEN; i++){
        for(j = 0; j < LEN; j++){
            matriz[i][j] = i + j;
        }
    }
}

void clean(int matriz[LEN][LEN]){
    int i, j;
    for(i = 0; i < LEN; i++){
        for(j = 0; j < LEN; j++){
            matriz[i][j] = 0;
        }
    }
}

void multiplica_matriz(int matrizA[LEN][LEN], int matrizB[LEN][LEN], int matrizC[LEN][LEN]){
    clean(matrizC);
    int i, j, k, r;
    for(i = 0; i < LEN; i++){
        for (k = 0; k < LEN; k++){
            r = matrizA[i][k];
            for (j = 0; j < LEN; j++){
                matrizC[i][j] += r * matrizB[k][j];
            }
        }
    }
}

void multiplica_matriz_slow(int matrizA[LEN][LEN], int matrizB[LEN][LEN], int matrizC[LEN][LEN]){
    clean(matrizC);
    int i, j, k, r;
    for(j = 0; j < LEN; j++){
        for (k = 0; k < LEN; k++){
            r = matrizB[k][j];
            for (i = 0; i < LEN; i++){
                matrizC[i][j] += r * matrizA[i][k];
            }
        }
    }
}

void print_matriz(int matriz[LEN][LEN]){
    int i, j;
    printf("{\n");
    for(i = 0; i < LEN; i++){
        printf("{");
        for(j = 0; j < LEN; j++){
            printf("%d", matriz[i][j]);
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
