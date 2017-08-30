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

static m *plus(m *restrict m_A, m *restrict m_B) {
    int n = m_A->re - m_A->rs + 1;
    m *restrict res = (m *)(malloc(sizeof(m)));
    res->rs = res->cs = 0;
    res->re = res->ce = n - 1;
    res->d = (double *)(malloc(n * n * sizeof(double)));

    int A_i = m_A->rs;
    int A_j = m_A->cs;
    int B_i = m_B->rs;
    int B_j = m_B->cs;

    double *restrict A = m_A->d;
    double *restrict B = m_B->d;
    double *restrict C = res->d;

    for (int i = 0; A_i <= m_A->re; A_i++, B_i++, i++) {
        for (int j = 0; A_j <= m_A->ce; A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] + B[B_i * n + B_j];
        }
    }

    return res;
}

static m *minus(m *restrict m_A, m *restrict m_B) {
    int n = m_A->re - m_A->rs + 1;
    m *restrict res = (m *)(malloc(sizeof(m)));
    res->rs = res->cs = 0;
    res->re = res->ce = n - 1;
    res->d = (double *)(malloc(n * n * sizeof(double)));

    int A_i = m_A->rs;
    int A_j = m_A->cs;
    int B_i = m_B->rs;
    int B_j = m_B->cs;

    double *restrict A = m_A->d;
    double *restrict B = m_B->d;
    double *restrict C = res->d;

    for (int i = 0; A_i <= m_A->re; A_i++, B_i++, i++) {
        for (int j = 0; A_j <= m_A->ce; A_j++, B_j++, j++) {
            C[i * n + j] = A[A_i * n + A_j] - B[B_i * n + B_j];
        }
    }

    return res;
}

static m *multiply(m *restrict m_A, m *restrict m_B) {
    int n = m_A->re - m_A->rs + 1;
    if (n <= 2) {
        double a, b, c, d, e, f, g, h;
        m *restrict m_C = (m *)(malloc(sizeof(m)));
        m_C->rs = m_A->rs;
        m_C->re = m_A->re;
        m_C->cs = m_A->cs;
        m_C->ce = m_A->ce;
        m_C->d = (double *)(malloc(n * n * sizeof(double)));

        double *restrict A = m_A->d;
        double *restrict B = m_B->d;
        double *restrict C = m_C->d;

        int A_line_start_idx = m_A->rs * n;
        int A_line_end_idx = m_A->rs * n + n - 1;
        int B_line_start_idx = m_B->rs * n;
        int B_line_end_idx = m_B->rs * n + n - 1;
        int A_cs_idx = m_A->cs;
        int A_ce_idx = m_A->cs + 1;
        int B_cs_idx = m_B->cs;
        int B_ce_idx = m_B->cs + 1;

        a = A[A_line_start_idx + A_cs_idx];
        b = A[A_line_start_idx + A_ce_idx];
        c = A[A_line_end_idx + A_cs_idx];
        d = A[A_line_end_idx + A_ce_idx];
        e = B[B_line_start_idx + B_cs_idx];
        f = B[B_line_start_idx + B_ce_idx];
        g = B[B_line_end_idx + B_cs_idx];
        h = B[B_line_end_idx + B_ce_idx];

        C[A_line_start_idx + A_cs_idx] = a * e + b * g;
        C[A_line_start_idx + A_ce_idx] = a * f + b * h;
        C[A_line_end_idx + A_cs_idx] = c * e + d * g;
        C[A_line_end_idx + A_ce_idx] = c * f + d * h;

        return m_C;
    }

    double *restrict d_A = m_A->d;
    double *restrict d_B = m_B->d;

    m *restrict A, *restrict B, *restrict C, *restrict D, *restrict E,
        *restrict F, *restrict G, *restrict H;
    m *restrict P1, *restrict P2, *restrict P3, *restrict P4, *restrict P5,
        *restrict P6, *restrict P7;
    m *restrict Q1, *restrict Q2, *restrict Q3, *restrict Q4;
    m *restrict res;
    int A_i, A_j;
    int i, j;

    A = (m *)(malloc(sizeof(m)));
    B = (m *)(malloc(sizeof(m)));
    C = (m *)(malloc(sizeof(m)));
    D = (m *)(malloc(sizeof(m)));
    E = (m *)(malloc(sizeof(m)));
    F = (m *)(malloc(sizeof(m)));
    G = (m *)(malloc(sizeof(m)));
    H = (m *)(malloc(sizeof(m)));
    res = (m *)(malloc(sizeof(m)));
    res->d = (double *)(malloc(n * n * sizeof(double)));

    A->d = B->d = C->d = D->d = d_A;
    E->d = F->d = G->d = H->d = d_B;

    res->rs = res->cs = 0;
    res->re = res->ce = n - 1;

    A->rs = m_A->rs;
    A->re = (m_A->re - m_A->rs) / 2 + m_A->rs;
    A->cs = m_A->cs;
    A->ce = (m_A->ce - m_A->cs) / 2 + m_A->cs;

    B->rs = m_A->rs;
    B->re = (m_A->re - m_A->rs) / 2 + m_A->rs;
    B->cs = (m_A->ce - m_A->cs) / 2 + m_A->cs + 1;
    B->ce = m_A->ce;

    C->rs = (m_A->re - m_A->rs) / 2 + m_A->rs + 1;
    C->re = m_A->re;
    C->cs = m_A->cs;
    C->ce = (m_A->ce - m_A->cs) / 2 + m_A->cs;

    D->rs = (m_A->re - m_A->rs) / 2 + m_A->rs + 1;
    D->re = m_A->re;
    D->cs = (m_A->ce - m_A->cs) / 2 + m_A->cs + 1;
    D->ce = m_A->ce;

    E->rs = m_B->rs;
    E->re = (m_B->re - m_B->rs) / 2 + m_B->rs;
    E->cs = m_B->cs;
    E->ce = (m_B->ce - m_B->cs) / 2 + m_B->cs;

    F->rs = m_B->rs;
    F->re = (m_B->re - m_B->rs) / 2 + m_B->rs;
    F->cs = (m_B->ce - m_B->cs) / 2 + m_B->cs + 1;
    F->ce = m_B->ce;

    G->rs = (m_B->re - m_B->rs) / 2 + m_B->rs + 1;
    G->re = m_B->re;
    G->cs = m_B->cs;
    G->ce = (m_B->ce - m_B->cs) / 2 + m_B->cs;

    H->rs = (m_B->re - m_B->rs) / 2 + m_B->rs + 1;
    H->re = m_B->re;
    H->cs = (m_B->ce - m_B->cs) / 2 + m_B->cs + 1;
    H->ce = m_B->ce;

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

    double *restrict m_C = res->d;

    int size = Q1->re - Q1->rs + 1;

    for (A_i = Q1->rs, i = 0; A_i <= Q1->re; A_i++, i++) {
        for (A_j = Q1->cs, j = 0; A_j <= Q1->ce; A_j++, j++) {
            m_C[A_i * n + A_j] = Q1->d[i * size + j];
        }
    }

    for (A_i = Q2->rs, i = 0; A_i <= Q2->re; A_i++, i++) {
        for (A_j = Q2->cs, j = n / 2; A_j <= Q2->ce; A_j++, j++) {
            m_C[A_i * n + A_j] = Q2->d[i * size + j];
        }
    }

    for (A_i = Q3->rs, i = n / 2; A_i <= Q3->re; A_i++, i++) {
        for (A_j = Q3->cs, j = 0; A_j <= Q3->ce; A_j++, j++) {
            m_C[i * n + j] = Q3->d[i * size + j];
        }
    }

    for (A_i = Q4->rs, i = n / 2; A_i <= Q4->re; A_i++, i++) {
        for (A_j = Q4->cs, j = n / 2; A_j <= Q4->ce; A_j++, j++) {
            m_C[A_i * n + A_j] = Q4->d[i * size + j];
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

    return res;
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm(int n, double *restrict mat_A, double *restrict mat_B,
                  double *restrict mat_C) {
    int lda = next_power_of_two(n);
    double *restrict d_A = (double *)(malloc(lda * lda * sizeof(double)));
    double *restrict d_B = (double *)(malloc(lda * lda * sizeof(double)));
    for (int i = 0; i < lda; i++) {
        for (int j = 0; j < lda; j++) {
            int pos_lda = i * lda + j;
            if (i < n && j < n) {
                int pos_n = i * n + j;
                d_A[pos_lda] = mat_A[pos_n];
                d_B[pos_lda] = mat_B[pos_n];
            } else {
                d_A[pos_lda] = 0.0;
                d_B[pos_lda] = 0.0;
            }
        }
    }
    m m_A = {0, lda - 1, 0, lda - 1, d_A};
    m m_B = {0, lda - 1, 0, lda - 1, d_B};
    m *m_C = multiply(&m_A, &m_B);

    // free(m_A.d);
    // free(m_B.d);

    double *restrict d_C = m_C->d;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int pos_n = i * n + j;
            int pos_lda = i * lda + j;
            mat_C[pos_n] = d_C[pos_lda];
        }
    }

    // free(m_C->d);
    // free(m_C);
}
