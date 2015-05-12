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

static double_complex woods_saxon(double_complex R,
                                  const struct woods_saxon_params *p)
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
        sf_bessel_jl((int)p->l, p->k1 * R) *
        sf_bessel_jl((int)p->l, p->k2 * R) *
        woods_saxon(R, p->ws);
}

double dmax(double x, double y)
{
    return x > y ? x : y;
}

void generate_v_matrix(double_complex *v, const double_complex *k,
                       size_t k_len, unsigned l,
                       const struct woods_saxon_params ws,
                       double r_max, double abs_err, double rel_err,
                       unsigned limit)
{
    unsigned total_num_eval = 0;
    for (size_t i = 0; i != k_len; ++i)
    for (size_t j = i; j != k_len; ++j) {
        const struct integrand_params p = {l, k[i], k[j], &ws};
        double est_err;
        unsigned num_eval;
        const double_complex integral =
            gk_cquad(integrand, &p, 0, r_max, abs_err, rel_err,
                     limit, &est_err, &num_eval);
        const double_complex element = 2 / M_PI * integral;
        v[i * k_len + j] = element;
        v[j * k_len + i] = element;
        if (!(est_err < dmax(abs_err, rel_err * cabs(integral)))) {
            fprintf(stderr, "WARNING: integral did not converge: "
                    "%i evals; error is %g, required %g\n",
                    num_eval, est_err, abs_err);
            fflush(stderr);
        }
        total_num_eval += num_eval;
    }
    fprintf(stderr, "note: number of interaction evaluations: %i evals\n",
            total_num_eval);
    fflush(stderr);
}
