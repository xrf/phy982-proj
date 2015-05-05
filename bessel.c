#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "bessel.h"

void zbesj_(const double *zr, const double *zi, const double *fnu,
            const int *kode, const int *m,
            double *cyr, double *cyi, int *nz, int *ierr);

int bessel_jl_e(double_complex *r, double m, double_complex z)
{
    static const int kode = 1, count = 1;
    const double zr = creal(z), zi = cimag(z);
    double rr, ri;
    int nz, ierr;
    zbesj_(&zr, &zi, &m, &kode, &count, &rr, &ri, &nz, &ierr);
    *r = rr + ri * _Complex_I;
    return ierr;
}

int sf_bessel_jl_e(double_complex *r, double l, double_complex z)
{
    const int e = bessel_jl_e(r, l + 0.5, z);
    switch (e) {
    case 0:
    case 3:
        *r = csqrt(M_PI_2 / z) * *r;
    }
    return e;
}

double_complex sf_bessel_jl(double l, double_complex z)
{
    double_complex r;
    const int e = sf_bessel_jl_e(&r, l, z);
    if (e) {
        fprintf(stderr, "sf_bessel_jl: error %d\n", e);
        abort();
    }
    return r;
}
