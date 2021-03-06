const char *dgemm_desc = "Simple blocked dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */

static void do_block(int lda, int M, int N, int K, double *restrict A, double *restrict B, double *restrict C) {
    /* For each row i of A */
    for (int i = 0; i < M; ++i) {
        /* For each column j of B */
        for (int j = 0; j < N; ++j) {
            /* Compute C(i,j) */
            int j_index = j * lda;
            int i_index = i + j_index;
            double cij = C[i_index];
            for (int k = 0; k < K; ++k) {
                double a = A[i + k * lda];
                double b = B[k + j_index];
                cij += a * b;
            }
            C[i_index] = cij;
        }
    }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int lda, double *restrict A, double *restrict B,
                  double *restrict C) {
    /* For each block-row of A */
    for (int i = 0; i < lda; i += BLOCK_SIZE) {
        /* For each block-column of B */
        for (int j = 0; j < lda; j += BLOCK_SIZE) {
            /* Accumulate block dgemms into block of C */
            for (int k = 0; k < lda; k += BLOCK_SIZE) {
                /* Correct block dimensions if block "goes off edge of" the
                 * matrix */
                int M = min(BLOCK_SIZE, lda - i);
                int N = min(BLOCK_SIZE, lda - j);
                int K = min(BLOCK_SIZE, lda - k);

                /* Perform individual block dgemm */
                do_block(lda, M, N, K, A + i + k * lda, B + k + j * lda,
                         C + i + j * lda);
            }
        }
    }
}
