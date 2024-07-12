#include <bessel.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef COMPARE_TO_FLINT
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <unistd.h>

typedef struct {
    acb_t res;
    acb_t z, nu;
    FILE* f;
    int fd;
} bessel_I_acb_handle_t;

static inline bessel_I_acb_handle_t* bessel_I_acb_handle_init( void ) {
    bessel_I_acb_handle_t* h = calloc(1, sizeof(bessel_I_acb_handle_t));
    acb_init(h->res);
    acb_init(h->z);
    acb_init(h->nu);
    h->f = tmpfile();
    h->fd = fileno(h->f);
    ftruncate(h->fd, 128);
    return h;
}

static inline void bessel_I_acb_handle_free( bessel_I_acb_handle_t* h) {
    acb_clear(h->res);
    acb_clear(h->z);
    acb_clear(h->nu);
    fclose(h->f);
    free(h);
}

static inline complex128_t bessel_I_acb( bessel_I_acb_handle_t* h, complex128_t nu, complex128_t z ) {
    complex128_t Cres;

    acb_set_d_d( h->nu, creal(nu), cimag(nu) );
    acb_set_d_d( h->z, creal(z), cimag(z) );
    acb_hypgeom_bessel_i( h->res, h->nu, h->z, 48 );

    rewind( h->f );
    arb_fprintn( h->f, acb_realref(h->res), 16, ARB_STR_NO_RADIUS );
    fputc(' ', h->f);
    arb_fprintn( h->f, acb_imagref(h->res), 16, ARB_STR_NO_RADIUS );

    rewind( h->f );
    fscanf( h->f, "%le %le", &Cres, (double*)(&Cres)+1);

    return Cres;
}

static inline complex128_t bessel_I_ratio_flint( bessel_I_acb_handle_t* h, complex128_t nu, complex128_t z ) {
    return bessel_I_acb(h, nu, z)/bessel_I_acb(h, nu-1.0, z);
}
#endif // COMPARE_TO_FLINT

int main() {
#ifdef COMPARE_TO_FLINT
    bessel_I_acb_handle_t* h = bessel_I_acb_handle_init();
#endif
    for (double g=-1.0; g<=1.0+1.e-7; g+=1.e-2)
    for (double c=-1.0; c<=1.0+1.e-7; c+=1.e-2) {
        complex128_t nu = 1.0 + I * g,
                     z = I * c;
        complex128_t J2 = bessel_I_ratio( nu, z, NULL ),
                     J3 = bessel_I_ratio_iter( nu, z, 100 );
#ifdef COMPARE_TO_FLINT
        complex128_t J1 = bessel_I_ratio_flint( h, nu, z );
        printf( "%.5e %.5e %.5e %.5e %.5e %.5e\n", creal(J1), cimag(J1), creal(J2), cimag(J2), creal(J3), cimag(J3) );
#else // COMPARE_TO_FLINT
        printf( "%.5e %.5e %.5e %.5e\n", creal(J2), cimag(J2), creal(J3), cimag(J3) );
#endif // COMPARE_TO_FLINT
    }
#ifdef COMPARE_TO_FLINT
    bessel_I_acb_handle_free(h);
#endif
}
