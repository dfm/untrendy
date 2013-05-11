#include <stdio.h>
#include <stdlib.h>

void fpcurf_(int *iopt, float x[], float y[], float w[], int *m, float *xb,
             float *xe, int *k, float *s, int *nest, float *tol, int *maxit,
             int *k1, int *k2, int *n, float t[], float c[], float *fp,
             float fpint[], float z[], float *a[], float *b[], float *g[],
             float *q[], int nrdata[], int *ier);


int max (int i1, int i2) {
    if (i1 >= i2) return i1;
    return i2;
}


int main()
{
    int i;
    int iopt = -1, m = 50, k = 3, k1 = k + 1, k2 = k + 2, nest = 15,
        n = nest, maxit = 20, ier;
    float s = -1.0, tol = 0.001, xb = -0.5, xe = 0.5, fp;

    float *x = malloc(m * sizeof(float)),
          *y = malloc(m * sizeof(float)),
          *w = malloc(m * sizeof(float)),
          *t = malloc(nest * sizeof(float)),
          *c = malloc(nest * sizeof(float)),
          *fpint = malloc(nest * sizeof(float)),
          *z = malloc(nest * sizeof(float)),
          **a = malloc(k1 * sizeof(float*)),
          **b = malloc(k2 * sizeof(float*)),
          **g = malloc(k2 * sizeof(float*)),
          **q = malloc(k1 * sizeof(float*));
    int *nrdata = malloc(nest * sizeof(int));

    for (i = 0; i < k1; ++i) {
        a[i] = malloc(nest * sizeof(float));
        g[i] = malloc(nest * sizeof(float));
    }

    for (i = 0; i < k2; ++i) {
        b[i] = malloc(nest * sizeof(float));
        q[i] = malloc(m * sizeof(float));
    }

    for (i = 0; i < m; ++i) {
        x[i] = xb + (xe - xb) * (float)i / (float)(m - 1);
        y[i] = 0.1 * x[i] * x[i] - 5;
        w[i] = 1.0;
    }

    for (i = k1; i < nest - k1; ++i) {
        t[i] = xb + (xe - xb) * (float)(i - k) / (float)(nest - k1);
        printf("%f\n", t[i]);
    }

    for (i = 0; i < k1; ++i) {
        t[i] = xb;
        t[nest - i - 1] = xe;
    }

    fpcurf_(&iopt, x, y, w, &m, &xb, &xe, &k, &s, &nest, &tol, &maxit,
            &k1, &k2, &n, t, c, &fp, fpint, z, a, b, g, q, nrdata, &ier);
    printf("%d\n", n);

    for (i = 0; i < nest - k1; ++i) {
        printf("%f %f\n", t[i + k - 1], c[i]);
        /* free(a[i]); */
        /* free(b[i]); */
        /* free(g[i]); */
    }

    /* for (i = 0; i < m; ++i) */
    /*     free(q[i]); */

    free(a);
    free(b);
    free(q);
    free(g);

    free(x);
    free(y);
    free(w);
    free(t);
    free(c);
    free(fpint);
    free(z);
    free(nrdata);

    return 0;
}
