#include <math.h>
#include <stdlib.h>

const char *dgemm_desc = "Strassen divide and conquer dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

typedef struct _matrix {
    int row_start, row_end, column_start, column_end;
    double *data;
} matrix;

static int next_power_of_two(int n) {
    return (int)(pow(2, ceil(log((double)n) / log(2.0))));
}

static matrix *plus(matrix *restrict matrix_A, matrix *restrict matrix_B) {
    matrix *restrict result = (matrix *)(malloc(sizeof(matrix)));
    int A_i = matrix_A->row_start;
    int A_j = matrix_A->column_start;
    int B_i = matrix_B->row_start;
    int B_j = matrix_B->column_start;
    int n = matrix_A->row_end - matrix_A->row_start + 1;
    result->data = (double *)(malloc(n * n * sizeof(double)));

    result->row_start = result->column_start = 0;
    result->row_end = result->column_end = n - 1;

    double *restrict A = matrix_A->data;
    double *restrict B = matrix_B->data;
    double *restrict C = result->data;

    for (int i = 0; A_i < matrix_A->row_end; A_i++, B_i++, i++) {
        for (int j = 0; A_j < matrix_A->column_end; A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] + B[B_i * n + B_j];
        }
    }

    return result;
}

static matrix *minus(matrix *restrict matrix_A, matrix *restrict matrix_B) {
    matrix *restrict result = (matrix *)(malloc(sizeof(matrix)));
    int A_i = matrix_A->row_start;
    int A_j = matrix_A->column_start;
    int B_i = matrix_B->row_start;
    int B_j = matrix_B->column_start;
    int n = matrix_A->row_end - matrix_A->row_start + 1;
    result->data = (double *)(malloc(n * n * sizeof(double)));

    result->row_start = result->column_start = 0;
    result->row_end = result->column_end = n - 1;

    double *restrict A = matrix_A->data;
    double *restrict B = matrix_B->data;
    double *restrict C = result->data;

    for (int i = 0; A_i < matrix_A->row_end; A_i++, B_i++, i++) {
        for (int j = 0; A_j < matrix_A->column_end; A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] - B[B_i * n + B_j];
        }
    }

    return result;
}

static matrix *multiply(matrix *restrict matrix_A, matrix *restrict matrix_B) {
    int n = matrix_A->row_end - matrix_A->row_start + 1;
    if (n == 2) {
        double a, b, c, d, e, f, g, h;
        matrix *matrix_C = (matrix *)(malloc(sizeof(matrix)));
        matrix_C->data = (double *)(malloc(sizeof(n * n * sizeof(double))));

        double *restrict A = matrix_A->data;
        double *restrict B = matrix_B->data;
        double *restrict C = matrix_C->data;

        int A_line_start_idx = matrix_A->row_start * n;
        int A_line_end_idx = A_line_start_idx + n;
        int B_line_start_idx = matrix_B->row_start * n;
        int B_line_end_idx = B_line_start_idx + n;
        int A_column_start_idx = matrix_A->column_start;
        int A_column_end_idx = matrix_A->column_start + 1;
        int B_column_start_idx = matrix_B->column_start;
        int B_column_end_idx = matrix_B->column_start + 1;

        a = A[A_line_start_idx + A_column_start_idx];
        b = A[A_line_start_idx + A_column_end_idx];
        c = A[A_line_end_idx + A_column_start_idx];
        d = A[A_line_end_idx + A_column_end_idx];
        e = B[B_line_start_idx + B_column_start_idx];
        f = B[B_line_start_idx + B_column_end_idx];
        g = B[B_line_end_idx + B_column_start_idx];
        h = B[B_line_end_idx + B_column_end_idx];

        C[A_line_start_idx + A_column_start_idx] = a * e + b * g;
        C[A_line_start_idx + A_column_end_idx] = a * f + b * h;
        C[A_line_end_idx + A_column_start_idx] = c * e + d * g;
        C[A_line_end_idx + A_column_end_idx] = c * f + d * h;

        return matrix_C;
    }

    double *restrict data_A = matrix_A->data;
    double *restrict data_B = matrix_B->data;

    matrix *restrict A, *restrict B, *restrict C, *restrict D, *restrict E,
        *restrict F, *restrict G, *restrict H;
    matrix *restrict P1, *restrict P2, *restrict P3, *restrict P4, *restrict P5,
        *restrict P6, *restrict P7;
    matrix *restrict Q1, *restrict Q2, *restrict Q3, *restrict Q4;
    matrix *restrict result;
    int A_i, A_j;
    int i, j;

    A = (matrix *)(malloc(sizeof(matrix)));
    B = (matrix *)(malloc(sizeof(matrix)));
    C = (matrix *)(malloc(sizeof(matrix)));
    D = (matrix *)(malloc(sizeof(matrix)));
    E = (matrix *)(malloc(sizeof(matrix)));
    F = (matrix *)(malloc(sizeof(matrix)));
    G = (matrix *)(malloc(sizeof(matrix)));
    H = (matrix *)(malloc(sizeof(matrix)));
    result = (matrix *)(malloc(sizeof(matrix)));
    result->data = (double *)(malloc(n * n * sizeof(double)));

    A->data = B->data = C->data = D->data = data_A;
    E->data = F->data = G->data = H->data = data_B;

    result->row_start = result->column_start = 0;
    result->row_end = result->column_end = n - 1;

    A->row_start = matrix_A->row_start;
    A->row_end = matrix_A->row_end / 2;
    A->column_start = matrix_A->column_start;
    A->column_end = matrix_A->column_end / 2;

    B->row_start = matrix_A->row_start;
    B->row_end = matrix_A->row_end / 2;
    B->column_start = matrix_A->column_end / 2 + 1;
    B->column_end = matrix_A->column_end;

    C->row_start = matrix_A->row_end / 2 + 1;
    C->row_end = matrix_A->row_end;
    C->column_start = matrix_A->column_start;
    C->column_end = matrix_A->column_end / 2;

    D->row_start = matrix_A->row_end / 2 + 1;
    D->row_end = matrix_A->row_end;
    D->column_start = matrix_A->column_end / 2 + 1;
    D->column_end = matrix_A->column_end;

    E->row_start = matrix_B->row_start;
    E->row_end = matrix_B->row_end / 2;
    E->column_start = matrix_B->column_start;
    E->column_end = matrix_B->column_end / 2;

    F->row_start = matrix_B->row_start;
    F->row_end = matrix_B->row_end / 2;
    F->column_start = matrix_B->column_end / 2 + 1;
    F->column_end = matrix_B->column_end;

    G->row_start = matrix_B->row_end / 2 + 1;
    G->row_end = matrix_B->row_end;
    G->column_start = matrix_B->column_start;
    G->column_end = matrix_B->column_end / 2;

    H->row_start = matrix_B->row_end / 2 + 1;
    H->row_end = matrix_B->row_end;
    H->column_start = matrix_B->column_end / 2 + 1;
    H->column_end = matrix_B->column_end;

    P1 = multiply(A, minus(F, H));
    P2 = multiply(plus(A, B), H);
    P3 = multiply(plus(C, D), E);
    P4 = multiply(D, minus(G, E));
    P5 = multiply(plus(A, D), plus(E, H));
    P6 = multiply(minus(B, D), plus(G, H));
    P7 = multiply(minus(A, C), plus(E, F));

    Q1 = plus(minus(plus(P5, P4), P2), P6);
    Q2 = plus(P1, P2);
    Q3 = plus(P3, P4);
    Q4 = minus(minus(plus(P1, P5), P3), P7);

    double *restrict matriz_C = result->data;

    for (m1_i = Q1->row_start, i = 0; m1_i <= Q1->row_end; m1_i++, i++)
        for (m1_j = Q1->column_start, j = 0; m1_j <= Q1->column_end;
             m1_j++, j++)
            matriz_C[i * n + j] = Q1->data[m1_i][m1_j];

    for (m1_i = Q2->row_start, i = 0; m1_i <= Q2->row_end; m1_i++, i++)
        for (m1_j = Q2->column_start, j = n / 2; m1_j <= Q2->column_end;
             m1_j++, j++)
            matriz_C[i * n + j] = Q2->data[m1_i][m1_j];

    for (m1_i = Q3->row_start, i = n / 2; m1_i <= Q3->row_end; m1_i++, i++)
        for (m1_j = Q3->column_start, j = 0; m1_j <= Q3->column_end;
             m1_j++, j++)
            matriz_C[i * n + j] = Q3->data[m1_i][m1_j];

    for (m1_i = Q4->row_start, i = n / 2; m1_i <= Q4->row_end; m1_i++, i++)
        for (m1_j = Q4->column_start, j = n / 2; m1_j <= Q4->column_end;
             m1_j++, j++)
            matriz_C[i * n + j] = Q4->data[m1_i][m1_j];

    return result;
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, double *restrict mat_A, double *restrict mat_B,
                  double *restrict mat_C) {
    int lda = next_power_of_two(n);
    matrix matrix_A = {0, lda - 1, 0, lda - 1, mat_A};
    matrix matrix_B = {0, lda - 1, 0, lda - 1, mat_B};
}
