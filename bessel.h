#pragma once

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> complex128_t;
extern "C" {
#else
#include <complex.h>
typedef _Complex double complex128_t;
#endif

// J_nu(z)
complex128_t bessel_J_complex( complex128_t nu, complex128_t z, int* info_ );
// Y_nu(z)
complex128_t bessel_Y_complex( complex128_t nu, complex128_t z, int* info_ );
// H^(1)_nu(z)
complex128_t hankel_1_complex( complex128_t nu, complex128_t z, int* info_ );
// H^(2)_nu(z)
complex128_t hankel_2_complex( complex128_t nu, complex128_t z, int* info_ );
// I_nu(z)
complex128_t bessel_I_complex( complex128_t nu, complex128_t z, int* info_ );
// I_nu(z) / I_nu-1(z)
complex128_t bessel_I_ratio( complex128_t nu, complex128_t z, int* info_ );
// I_nu(z) / I_nu-1(z), but with continued fraction implementation of imax
// iterations. When imax <= 0, fall back to bessel_I_ratio()
complex128_t bessel_I_ratio_iter( complex128_t nu, complex128_t z, int imax );

#ifdef __cplusplus
}
#endif
