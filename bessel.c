#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bessel.h"

void zbesj_(const double *zr, const double *zi, const double *fnu,
            const int *kode, const int *m,
            double *cyr, double *cyi, int *nz, int *ierr);

int bessel_jl_e(double_complex *r, double m, double_complex z)
{
    static const int kode = 1, count = 1;
    const double zr = creal(z), zi = cimag(z);
    double rr = NAN, ri = NAN;
    int nz, ierr;
    zbesj_(&zr, &zi, &m, &kode, &count, &rr, &ri, &nz, &ierr);
    *r = rr + ri * _Complex_I;
    return ierr;
}

int spherical_bessel_jl_e(double_complex *r, double l, double_complex z)
{
    const int e = bessel_jl_e(r, l + .5, z);
    *r = csqrt(M_PI_2 / z) * *r;
    return e;
}

double_complex spherical_bessel_jl(double l, double_complex z)
{
    double_complex r;
    const int e = spherical_bessel_jl_e(&r, l, z);
    switch (e) {
    case 0:
        break;
    case 1:
        fprintf(stderr, "spherical_bessel_jl(%.17g, (%.17g)+(%.17g)i) "
                "[error]: input error\n", l, creal(z), cimag(z));
        abort();
    case 2:
        fprintf(stderr, "spherical_bessel_jl(%.17g, (%.17g)+(%.17g)i) "
                "[warning]: overflow, Im(z) too large\n",
                l, creal(z), cimag(z));
        break;
    case 3:
        fprintf(stderr, "spherical_bessel_jl(%.17g, (%.17g)+(%.17g)i) "
                "[note]: |z| or n large, losses of significance by "
                "argument reduction process less than half of machine "
                "accuracy\n", l, creal(z), cimag(z));
        break;
    case 4:
        fprintf(stderr, "spherical_bessel_jl(%.17g, (%.17g)+(%.17g)i) "
                "[warning]: |z| or n too large, complete losses of "
                "significance by argument reduction\n",
                l, creal(z), cimag(z));
        break;
    case 5:
        fprintf(stderr, "spherical_bessel_jl(%.17g, (%.17g)+(%.17g)i) "
                "[warning]: algorithm termination condition not met\n",
                l, creal(z), cimag(z));
        break;
    default:
        fprintf(stderr, "spherical_bessel_jl "
                "[error]: unknown error code of %d\n", e);
        abort();
    }
    return r;
}

void spherical_bessel_jl_many(double_complex *j, double l,
                              const double_complex *z, size_t count)
{
    const double_complex *z_end = z + count;
    for (; z != z_end; ++j, ++z)
        *j = spherical_bessel_jl(l, *z);
}
