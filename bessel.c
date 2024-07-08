#include "bessel.h"
#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

extern void _FORTRAN_bessel1( const complex128_t* znu, const complex128_t* zz, complex128_t* zans, int* info );
extern void _FORTRAN_bessel2( const complex128_t* znu, const complex128_t* zz, complex128_t* zans, int* info );
extern void _FORTRAN_hankel1( const complex128_t* znu, const complex128_t* zz, complex128_t* zans, int* info );
extern void _FORTRAN_hankel2( const complex128_t* znu, const complex128_t* zz, complex128_t* zans, int* info );

complex128_t bessel_J_complex( complex128_t nu, complex128_t z, int* info_ ) {
    complex128_t result = 0;
    int info = 0;
    _FORTRAN_bessel1( &nu, &z, &result, &info );
    if (info_) *info_ = info;
    return result;
}

complex128_t bessel_Y_complex( complex128_t nu, complex128_t z, int* info_ ) {
    complex128_t result = 0;
    int info = 0;
    _FORTRAN_bessel2( &nu, &z, &result, &info );
    if (info_) *info_ = info;
    return result;
}

complex128_t hankel_1_complex( complex128_t nu, complex128_t z, int* info_ ) {
    complex128_t result = 0;
    int info = 0;
    _FORTRAN_hankel1( &nu, &z, &result, &info );
    if (info_) *info_ = info;
    return result;
}

complex128_t hankel_2_complex( complex128_t nu, complex128_t z, int* info_ ) {
    complex128_t result = 0;
    int info = 0;
    _FORTRAN_hankel2( &nu, &z, &result, &info );
    if (info_) *info_ = info;
    return result;
}

complex128_t bessel_I_complex( complex128_t nu, complex128_t z, int* info_ ) {
    return cexp( nu * M_PI * I * 0.5 ) * bessel_J_complex( nu, z * cexp( -M_PI * I * 0.5 ), info_ );
}

complex128_t bessel_I_ratio( complex128_t nu, complex128_t z, int* info_ ) {
    int i1 = 0, i2 = 0;
    complex128_t bInu = bessel_I_complex(nu, z, &i1);
    complex128_t bInu1 = bessel_I_complex(nu-1.0, z, &i2);
    if (info_) *info_ = i1 > i2 ? i1 : i2;
    return bInu/bInu1;
}

complex128_t bessel_I_ratio_iter( complex128_t nu, complex128_t z, int imax ) {
    if (imax <= 0)
        return bessel_I_ratio( nu, z, NULL );
    complex128_t val = 0.0;
    for (int i=imax; i>=0; i--) {
        val = 1./(val + 2.*(nu+(double)i)/z);
    }
    return val;
}

