#pragma once

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> complex128_t;
extern "C" {
#else
#include <complex.h>
typedef _Complex double complex128_t;
#endif

/*
 * The info argument may be NULL, in which case it is ignored. Else, the
 * behavior of mod_zbes.f90 is inherited:
 *
 * info=0:  normal output; relative error of zans is less than
 *          MAX(epsilon1,1E3*epsilon0)
 * info=10: zans is not reliable because of one of the following reasons.
 *          (1) There is a possibility of an overflow.
 *          (2) There is a possibility of an underflow.
 *          (3) The precision of zans is not sufficient.
 *          (4) The output zans is indefinite theoretically.
 *              Function zbessel1(zunit,0) is indefinite for example.
 * info=20: out of range. The cases of ABS(re_zz)>ai_arg_m for example.
 */

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
