#include <assert.h>
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bessel.h"
#include "gauss_kronrod.h"

struct woods_saxon_params {
    double V0, Rws, aws;
};

static double_complex woods_saxon(double R, const struct woods_saxon_params *p)
{
    return p->V0 / (1 + exp((R - p->Rws) / p->aws));
}

struct integrand_params {
    unsigned l;
    double_complex k1, k2;
    const struct woods_saxon_params *ws;
};

static double_complex integrand(double R, const void *ctx)
{
    const struct integrand_params *p =
        (const struct integrand_params *)ctx;
    return R * R *
        spherical_bessel_jl((int)p->l, p->k1 * R) *
        spherical_bessel_jl((int)p->l, p->k2 * R) *
        woods_saxon(R, p->ws);
}

static double dmax(double x, double y)
{
    return x > y ? x : y;
}

static double_complex v_matrix_element(unsigned *total_num_eval,
                                       const struct integrand_params *params,
                                       double r_max, double abs_err,
                                       double rel_err, unsigned limit)
{
    double est_err;
    unsigned num_eval;
    const double_complex integral = gk_cquad(integrand, params, 0, r_max,
                                             abs_err, rel_err, limit,
                                             &est_err, &num_eval);
    const double_complex element = 2 / M_PI * integral;
    if (!(est_err < dmax(abs_err, rel_err * cabs(integral)))) {
        fprintf(stderr, "WARNING: integral did not converge: "
                "%i evals; error is %g, required %g\n",
                num_eval, est_err, abs_err);
        fflush(stderr);
    }
    if (total_num_eval)
        *total_num_eval += num_eval;
    return element;
}

/* note: k2 can be NULL (with k2_len ignored),
         in which case k2 is assumed to be same as k1 */
unsigned generate_v_matrix(double_complex *v, size_t v_stride,
                           const double_complex *k1, size_t k1_len,
                           const double_complex *k2, size_t k2_len,
                           unsigned l, struct woods_saxon_params ws,
                           double r_max, double abs_err, double rel_err,
                           unsigned limit)
{
    unsigned total_num_eval = 0;
    size_t i, j;
    assert(v);
    assert(k1);
    if (k2) {                           /* full matrix */
        for (i = 0; i != k1_len; ++i)
        for (j = 0; j != k2_len; ++j) {
            const struct integrand_params p = {l, k1[i], k2[j], &ws};
            v[i * v_stride + j] =
                v_matrix_element(&total_num_eval, &p, r_max,
                                 abs_err, rel_err, limit);
        }
    } else {                            /* upper triangular part */
        for (i = 0; i != k1_len; ++i)
        for (j = i; j != k1_len; ++j) {
            const struct integrand_params p = {l, k1[i], k1[j], &ws};
            const double_complex element =
                v_matrix_element(&total_num_eval, &p, r_max,
                                 abs_err, rel_err, limit);
            v[i * v_stride + j] = element;
            v[j * v_stride + i] = element;
        }
    }
    return total_num_eval;
}
