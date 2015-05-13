#ifndef G_JJVGRG0NJNI85L90IRHNAWQSB5FCL
#define G_JJVGRG0NJNI85L90IRHNAWQSB5FCL
#include "math_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

int bessel_jl_e(double_complex *r, double m, double_complex z);

int spherical_bessel_jl_e(double_complex *r, double l, double_complex z);

double_complex spherical_bessel_jl(double l, double_complex z);

void spherical_bessel_jl_many(double_complex *j, double l,
                              const double_complex *z, size_t count);

#ifdef __cplusplus
}
#endif
#endif
