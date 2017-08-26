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

static inline double get(int n, int i, int j, matrix *restrict matrix_A) {
    return (i >= n || j >= n) ? 0.0 : matrix_A->data[i * n + j];
}

static matrix *plus(matrix *restrict matrix_A, matrix *restrict matrix_B) {
    matrix *result = (matrix *)(malloc(sizeof(matrix)));
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
    matrix *result = (matrix *)(malloc(sizeof(matrix)));
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
        matrix_c->data = (double *)(malloc(sizeof(n * n * sizeof(double))));

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

        a = A[a_index];
        b = A[A_line_start_idx + A_column_end_idx];
        c = A[A_line_end_idx + A_column_start_idx];
        d = A[A_line_end_idx + A_column_end_idx];
        e = B[B_line_start_idx + B_column_start_idx];
        f = B[B_line_start_idx + B_column_end_idx];
        g = B[B_line_end_idx + B_column_start_idx];
        h = B[B_line_end_idx + B_column_end_idx];

        C[a_index] = a * e + b * g;
        C[A_line_start_idx + A_column_end_idx] = a * f + b * h;
        C[A_line_end_idx + A_column_start_idx] = c * e + d * g;
        C[A_line_end_idx + A_column_end_idx] = c * f + d * h;

        return matrix_C;
    }

    matrix *A, *B, *C, *D, *E, *F, *G, *H;
    matrix *P1, *P2, *P3, *P4, *P5, *P6, *P7;
    matrix *Q1, *Q2, *Q3, *Q4;
    matrix *result;
    int A_i, A_j;
    int i, j;

    result.row_start = result.column_start = 0;
    result.row_end = result.column_end = n - 1;
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
