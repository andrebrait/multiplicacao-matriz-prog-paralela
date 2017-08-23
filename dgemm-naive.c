const char *dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, double *restrict A, double *restrict B,
                  double *restrict C) {
    /* For each row i of A */
    for (int i = 0; i < n; ++i) {
        /* For each column j of B */
        for (int j = 0; j < n; ++j) {
            /* Compute C(i,j) */
            int j_index = j * n;
            int i_index = i + j_index;
            double cij = C[i_index];
            for (int k = 0; k < n; k++) {
                double a = A[i + k * n];
                double b = B[k + j_index];
                cij += a * b;
            }
            C[i_index] = cij;
        }
    }
}
