#ifndef G_JW36MOID287QXQYBHSFCVY6VV6DCI
#define G_JW36MOID287QXQYBHSFCVY6VV6DCI
#include "math_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

double gk_quad(double (*f)(double x, const void *ctx),
               const void *ctx, double a, double b,
               double abs_err, double rel_err, unsigned limit,
               double *est_err, unsigned *num_eval);

double_complex gk_cquad(double_complex (*f)(double x, const void *ctx),
                       const void *ctx, double a, double b,
                       double abs_err, double rel_err, unsigned limit,
                       double *est_err, unsigned *num_eval);
#ifdef __cplusplus
}
#endif
#endif
