#include <stdlib.h>
#include <math.h>
#include "untrendy.h"

double kernel (int n, double *t, double t0, double dt, double *softr)
{
    int i;
    double norm = 0.0, val = 0.0, k, delta;

    for (i = 0; i < n; ++i) {
        delta = t[i] - t0;
        if (delta >= -dt && delta <= dt) {

            if (delta >= 0.0) {
                k = delta / dt - 1.0;
                k *= k;
            } else {
                k = delta / dt + 1.0;
                k *= -k;
            }

            norm += k * k;
            val += k * softr[i];
        }
    }
    val /= norm;
    return val * val;
}

int find_discontinuities (int n, double *t, double *chi, double dt, double Q,
                          double thresh)
{
    int i, ind = -1;
    double tmid, val, maxv = 0.0;
    double *softr = (double*)malloc(n * sizeof(double));

    for (i = 0; i < n; ++i)
        softr[i] = sqrt(chi[i] * Q / (Q + chi[i] * chi[i]));

    for (i = 0; i < n - 1; ++i) {
        tmid = 0.5 * (t[i] + t[i + 1]);
        val = kernel(n, t, tmid, dt, softr);
        if (val >= thresh && val >= maxv) {
            ind = i;
            maxv = val;
        }
    }

    free(softr);
    return ind;
}
