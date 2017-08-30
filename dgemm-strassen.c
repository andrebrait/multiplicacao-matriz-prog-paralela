#include <math.h>
#include <stdlib.h>

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

const char *dgemm_desc = "Strassen divide and conquer dgemm.";

typedef struct _m {
    int rs, re, cs, ce;
    double *d;
} m;

static int next_power_of_two(int n) {
    return (int)(pow(2, ceil(log((double)n) / log(2.0))));
}

static m *plus(m *restrict mat_A, m *restrict mat_B) {
    int i, j, A_i, A_j, B_i, B_j;

    int n = mat_A->re - mat_A->rs + 1;
    m *restrict mat_C = (m *)(malloc(sizeof(m)));
    mat_C->rs = mat_C->cs = 0;
    mat_C->re = mat_C->ce = n - 1;
    mat_C->d = (double *)(malloc(n * n * sizeof(double)));

    double *restrict A = mat_A->d;
    double *restrict B = mat_B->d;
    double *restrict C = mat_C->d;

    for (i = 0, A_i = mat_A->rs, B_i = mat_B->rs; A_i <= mat_A->re;
         A_i++, B_i++, i++) {
        for (j = 0, A_j = mat_A->cs, B_j = mat_B->cs; A_j <= mat_A->ce;
             A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] + B[B_i * n + B_j];
        }
    }

    return mat_C;
}

static m *minus(m *restrict mat_A, m *restrict mat_B) {
    int i, j, A_i, A_j, B_i, B_j;

    int n = mat_A->re - mat_A->rs + 1;
    m *restrict mat_C = (m *)(malloc(sizeof(m)));
    mat_C->rs = mat_C->cs = 0;
    mat_C->re = mat_C->ce = n - 1;
    mat_C->d = (double *)(malloc(n * n * sizeof(double)));

    double *restrict A = mat_A->d;
    double *restrict B = mat_B->d;
    double *restrict C = mat_C->d;

    for (i = 0, A_i = mat_A->rs, B_i = mat_B->rs; A_i <= mat_A->re;
         A_i++, B_i++, i++) {
        for (j = 0, A_j = mat_A->cs, B_j = mat_B->cs; A_j <= mat_A->ce;
             A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] - B[B_i * n + B_j];
        }
    }

    return mat_C;
}

static m *multiply(m *restrict mat_A, m *restrict mat_B) {
    int n = mat_A->re - mat_A->rs + 1;
    m *mat_C = (m *)(malloc(sizeof(m)));
    mat_C->rs = mat_C->cs = 0;
    mat_C->re = mat_C->ce = n - 1;
    mat_C->d = (double *)(malloc(n * n * sizeof(double)));
    double *restrict m_A = mat_A->d;
    double *restrict m_B = mat_B->d;
    double *restrict m_C = mat_C->d;

    if (n == 2) {
        double a, b, c, d, e, f, g, h;

        a = m_A[mat_A->rs * n + mat_A->cs];
        b = m_A[mat_A->rs * n + mat_A->ce];
        c = m_A[mat_A->re * n + mat_A->cs];
        d = m_A[mat_A->re * n + mat_A->ce];
        e = m_B[mat_B->rs * n + mat_B->cs];
        f = m_B[mat_B->rs * n + mat_B->ce];
        g = m_B[mat_B->re * n + mat_B->cs];
        h = m_B[mat_B->re * n + mat_B->ce];

        m_C[0] = a * e + b * g;
        m_C[1] = a * f + b * h;
        m_C[2] = c * e + d * g;
        m_C[3] = c * f + d * h;
    } else if (n == 1) {
        m_C[0] = m_A[0] * m_B[0];
    } else {
        m *restrict A, *restrict B, *restrict C, *restrict D, *restrict E,
            *restrict F, *restrict G, *restrict H;
        m *restrict P1, *restrict P2, *restrict P3, *restrict P4, *restrict P5,
            *restrict P6, *restrict P7;
        m *restrict Q1, *restrict Q2, *restrict Q3, *restrict Q4;
        m *restrict res;
        int Q_i, Q_j;
        int i, j;

        A = (m *)(malloc(sizeof(m)));
        B = (m *)(malloc(sizeof(m)));
        C = (m *)(malloc(sizeof(m)));
        D = (m *)(malloc(sizeof(m)));
        E = (m *)(malloc(sizeof(m)));
        F = (m *)(malloc(sizeof(m)));
        G = (m *)(malloc(sizeof(m)));
        H = (m *)(malloc(sizeof(m)));

        A->d = B->d = C->d = D->d = m_A;
        E->d = F->d = G->d = H->d = m_B;

        A->rs = mat_A->rs;
        A->re = (mat_A->re - mat_A->rs) / 2 + mat_A->rs;
        A->cs = mat_A->cs;
        A->ce = (mat_A->ce - mat_A->cs) / 2 + mat_A->cs;

        B->rs = mat_A->rs;
        B->re = (mat_A->re - mat_A->rs) / 2 + mat_A->rs;
        B->cs = (mat_A->ce - mat_A->cs) / 2 + mat_A->cs + 1;
        B->ce = mat_A->ce;

        C->rs = (mat_A->re - mat_A->rs) / 2 + mat_A->rs + 1;
        C->re = mat_A->re;
        C->cs = mat_A->cs;
        C->ce = (mat_A->ce - mat_A->cs) / 2 + mat_A->cs;

        D->rs = (mat_A->re - mat_A->rs) / 2 + mat_A->rs + 1;
        D->re = mat_A->re;
        D->cs = (mat_A->ce - mat_A->cs) / 2 + mat_A->cs + 1;
        D->ce = mat_A->ce;

        E->rs = mat_B->rs;
        E->re = (mat_B->re - mat_B->rs) / 2 + mat_B->rs;
        E->cs = mat_B->cs;
        E->ce = (mat_B->ce - mat_B->cs) / 2 + mat_B->cs;

        F->rs = mat_B->rs;
        F->re = (mat_B->re - mat_B->rs) / 2 + mat_B->rs;
        F->cs = (mat_B->ce - mat_B->cs) / 2 + mat_B->cs + 1;
        F->ce = mat_B->ce;

        G->rs = (mat_B->re - mat_B->rs) / 2 + mat_B->rs + 1;
        G->re = mat_B->re;
        G->cs = mat_B->cs;
        G->ce = (mat_B->ce - mat_B->cs) / 2 + mat_B->cs;

        H->rs = (mat_B->re - mat_B->rs) / 2 + mat_B->rs + 1;
        H->re = mat_B->re;
        H->cs = (mat_B->ce - mat_B->cs) / 2 + mat_B->cs + 1;
        H->ce = mat_B->ce;

        P1 = multiply(A, minus(F, H));
        P2 = multiply(plus(A, B), H);
        P3 = multiply(plus(C, D), E);
        P4 = multiply(D, minus(G, E));
        P5 = multiply(plus(A, D), plus(E, H));
        P6 = multiply(minus(B, D), plus(G, H));
        P7 = multiply(minus(A, C), plus(E, F));

        // free(A);
        // free(B);
        // free(C);
        // free(D);
        // free(E);
        // free(F);
        // free(G);

        Q1 = plus(minus(plus(P5, P4), P2), P6);
        Q2 = plus(P1, P2);
        Q3 = plus(P3, P4);
        Q4 = minus(minus(plus(P1, P5), P3), P7);

        // free(P1->d);
        // free(P2->d);
        // free(P3->d);
        // free(P4->d);
        // free(P5->d);
        // free(P6->d);
        // free(P7->d);
        // free(P1);
        // free(P2);
        // free(P3);
        // free(P4);
        // free(P5);
        // free(P6);
        // free(P7);

        int size = Q1->re - Q1->rs + 1;

        for (Q_i = Q1->rs, i = 0; Q_i <= Q1->re; Q_i++, i++) {
            for (Q_j = Q1->cs, j = 0; Q_j <= Q1->ce; Q_j++, j++) {
                m_C[i * n + j] = Q1->d[Q_i * size + Q_j];
            }
        }

        for (Q_i = Q2->rs, i = 0; Q_i <= Q2->re; Q_i++, i++) {
            for (Q_j = Q2->cs, j = n / 2; Q_j <= Q2->ce; Q_j++, j++) {
                m_C[i * n + j] = Q2->d[Q_i * size + Q_j];
            }
        }

        for (Q_i = Q3->rs, i = n / 2; Q_i <= Q3->re; Q_i++, i++) {
            for (Q_j = Q3->cs, j = 0; Q_j <= Q3->ce; Q_j++, j++) {
                m_C[i * n + j] = Q3->d[Q_i * size + Q_j];
            }
        }

        for (Q_i = Q4->rs, i = n / 2; Q_i <= Q4->re; Q_i++, i++) {
            for (Q_j = Q4->cs, j = n / 2; Q_j <= Q4->ce; Q_j++, j++) {
                m_C[i * n + j] = Q4->d[Q_i * size + Q_j];
            }
        }

        // free(Q1->d);
        // free(Q2->d);
        // free(Q3->d);
        // free(Q4->d);
        // free(Q1);
        // free(Q2);
        // free(Q3);
        // free(Q4);
    }
    return mat_C;
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, double *restrict A, double *restrict B,
                  double *restrict C) {
    int lda = next_power_of_two(n);
    double *restrict d_A = (double *)(malloc(lda * lda * sizeof(double)));
    double *restrict d_B = (double *)(malloc(lda * lda * sizeof(double)));
    for (int i = 0; i < lda; i++) {
        for (int j = 0; j < lda; j++) {
            int pos_lda = i * lda + j;
            if (i < n && j < n) {
                int pos_n = i * n + j;
                d_A[pos_lda] = A[pos_n];
                d_B[pos_lda] = B[pos_n];
            } else {
                d_A[pos_lda] = 0.0;
                d_B[pos_lda] = 0.0;
            }
        }
    }
    m mat_A = {0, lda - 1, 0, lda - 1, d_A};
    m mat_B = {0, lda - 1, 0, lda - 1, d_B};
    m *mat_C = multiply(&mat_A, &mat_B);

    // free(mat_A.d);
    // free(mat_B.d);

    double *restrict d_C = mat_C->d;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int pos_n = i * n + j;
            int pos_lda = i * lda + j;
            C[pos_n] = d_C[pos_lda];
        }
    }

    // free(m_C->d);
    // free(m_C);
}
