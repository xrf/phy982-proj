#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "xmalloc.h"
#include "gauss_kronrod.h"

void dqage_(double (*f)(const double *x), const double *a, const double *b,
            const double *epsabs, const double *epsrel, const int *key,
            const int *limit, double *result, double *abserr, int *neval,
            int *ier, double *alist, double *blist, double *rlist,
            double *elist, int *iord, int *last);

static double (*f_f)(double x, const void *ctx);
static const void *f_ctx;
static double f_fun(const double *x)
{
    return f_f(*x, f_ctx);
}

// doesn't seem to make much difference for our integrals
static int gk_key = 3;

// not thread-safe!
//  7 - 15 points if key < 2,
// 10 - 21 points if key = 2,
// 15 - 31 points if key = 3,
// 20 - 41 points if key = 4,
// 25 - 51 points if key = 5,
// 30 - 61 points if key > 5.
/* error <= max(epsabs, epsrel * abs(exact_integral)) */
double gk_quad(double (*f)(double x, const void *ctx),
               const void *ctx, double a, double b,
               double abs_err, double rel_err, unsigned limit,
               double *est_err, unsigned *num_eval)
{
    const int limit_ = (int)limit; // potential overflow
    double result = NAN, est_err_ = NAN;
    int num_eval_ = 0, ier = -1, last;
    double * work = (double *)xmalloc(sizeof(*work)  * limit * 4);
    int    *iwork = (int    *)xmalloc(sizeof(*iwork) * limit);
    f_f   = f;
    f_ctx = ctx;
    dqage_(f_fun, &a, &b, &abs_err, &rel_err, &gk_key, &limit_,
           &result, &est_err_, &num_eval_, &ier,
           work, work + limit, work + limit * 2, work + limit * 3,
           iwork, &last);
    if (ier) {
        fprintf(stderr, "gk_quad: [WARNING] ier = %i\n", ier);
        fflush(stderr);
    }
    if (est_err)
        *est_err = est_err_;
    if (num_eval)
        *num_eval = (unsigned)num_eval_;
    free(work);
    free(iwork);
    return result;
}

/****************************************************************************/

#undef SQRT_2
#define SQRT_2 1.414213562373095048801688724209698078569671875376948073176679

struct cf_ctx {
    double_complex (*f)(double x, const void *ctx);
    const void *ctx;
};

static double cf_real(double x, const void *ctx)
{
    const struct cf_ctx *p = (const struct cf_ctx *)ctx;
    return creal((*p->f)(x, p->ctx));
}

static double cf_imag(double x, const void *ctx)
{
    const struct cf_ctx *p = (const struct cf_ctx *)ctx;
    return cimag((*p->f)(x, p->ctx));
}

double_complex gk_cquad(double_complex (*f)(double x, const void *ctx),
                        const void *ctx, double a, double b,
                        double abs_err, double rel_err, unsigned limit,
                        double *est_err, unsigned *num_eval)
{
    const struct cf_ctx p = {f, ctx};
    double est_err_real, est_err_imag;
    unsigned num_eval_real, num_eval_imag;
    abs_err /= SQRT_2;
    rel_err /= SQRT_2;
    const double real = gk_quad(cf_real, &p, a, b, abs_err, rel_err, limit,
                                est_err  ? &est_err_real  : NULL,
                                num_eval ? &num_eval_real : NULL);
    const double imag = gk_quad(cf_imag, &p, a, b, abs_err, rel_err, limit,
                                est_err  ? &est_err_imag  : NULL,
                                num_eval ? &num_eval_imag : NULL);
    if (num_eval)
        *num_eval = num_eval_real + num_eval_imag;
    if (est_err)
        *est_err = cabs(est_err_real + _Complex_I * est_err_imag);
    return real + _Complex_I * imag;
}
