#ifndef G_5I16AJKS8G45ACO37R2IF93JRGGRC
#define G_5I16AJKS8G45ACO37R2IF93JRGGRC

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> double_complex;
#else
#if __STDC_VERSION__ >= 199901L
#include <complex.h>
typedef double _Complex double_complex;
#elif !defined(HAVE_DOUBLE_COMPLEX)
typedef struct { double real, imag; } double_complex;
#define HAVE_DOUBLE_COMPLEX
#endif
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058
#endif

#ifndef M_PI_2
#define M_PI_2 1.5707963267948966192313216916397514420985846996875529
#endif

#endif
