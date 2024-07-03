#include <stdio.h>
#include <bessel.h>

int main() {
    for (double g=-1.0; g<=1.0+1.e-7; g+=1.e-2)
    for (double c=-1.0; c<=1.0+1.e-7; c+=1.e-2) {
        complex128_t J1 = bessel_I_ratio( 1.0 + I * g, I * c, NULL ),
                     J2 = bessel_I_ratio_iter( 1.0 + I * g, I * c, 100 );
        printf( "%.5e%+.5ej %.5e%+.5ej\n", creal(J1), cimag(J1), creal(J2), cimag(J2) );
    }
}
