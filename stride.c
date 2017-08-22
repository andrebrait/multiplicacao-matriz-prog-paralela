#include<stdio.h>

#define LEN 1024
#define STRIDE 1

int vetor[LEN][LEN];

int main(int argc, char** argv){
    int i;
    for(i = 0; i < LEN; i += STRIDE){
        vetor[i] = 0;
    }
}