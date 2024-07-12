#include <bessel.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef COMPARE_TO_FLINT
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

static inline complex128_t bessel_I_acb( complex128_t nu_, complex128_t z_ ) {
    acb_t res;
    acb_init(res);

    acb_t z, nu;
    acb_init(z),
    acb_init(nu);
    acb_set_d_d( nu, creal(nu_), cimag(nu_) );
    acb_set_d_d( z, creal(z_), cimag(z_) );
    acb_hypgeom_bessel_i(res, nu, z, 48);
    acb_clear(z),
    acb_clear(nu);

    arb_t res_r, res_i;
    arb_init( res_r );
    arb_init( res_i );
    acb_get_real( res_r, res );
    acb_get_imag( res_i, res );
    acb_clear(res);

    char* r = arb_get_str( res_r, 16, ARB_STR_NO_RADIUS );
    double rr = atof(r);
    flint_free(r);
    char* i = arb_get_str( res_i, 16, ARB_STR_NO_RADIUS );
    double ii = atof(i);
    flint_free(i);

    arb_clear( res_r );
    arb_clear( res_i );

    return rr + I*ii;
}

static inline complex128_t bessel_I_ratio_flint( complex128_t nu, complex128_t z ) {
    return bessel_I_acb(nu, z)/bessel_I_acb(nu-1.0, z);
}
#endif // COMPARE_TO_FLINT

int main() {
    for (double g=-1.0; g<=1.0+1.e-7; g+=1.e-2)
    for (double c=-1.0; c<=1.0+1.e-7; c+=1.e-2) {
        complex128_t nu = 1.0 + I * g,
                     z = I * c;
        complex128_t J2 = bessel_I_ratio( nu, z, NULL ),
                     J3 = bessel_I_ratio_iter( nu, z, 100 );
#ifdef COMPARE_TO_FLINT
        complex128_t J1 = bessel_I_ratio_flint( nu, z );
        printf( "%.5e %.5e %.5e %.5e %.5e %.5e\n", creal(J1), cimag(J1), creal(J2), cimag(J2), creal(J3), cimag(J3) );
#else // COMPARE_TO_FLINT
        printf( "%.5e %.5e %.5e %.5e\n", creal(J2), cimag(J2), creal(J3), cimag(J3) );
#endif // COMPARE_TO_FLINT
    }
}
