# Cylindrical Function of Complex Order and Complex Argument in C

This is a wrapper around the algorithms published in [this
article](https://doi.org/10.1145/1916461.1916471). It computes the functions
```math
J_\nu(z),~ Y_\nu(z),~H_\nu^{(1)}(z),~H_\nu^{(2)}(z)
```
for $`\nu\in\mathbb C`$ and $`z\in\mathbb C`$. Based on the implementation of
$`J_\nu(z)`$, this package further includes an implementation of $`I_\nu(z)`$,
see [DLMF](https://dlmf.nist.gov/10.27). It also includes an implementation of
```math
\mathcal{J}_\nu(z) = \frac{I_\nu(z)}{I_{\nu-1}(z)}
```
in both the implementation via $`I_\nu(z)`$ as well as the [continued fraction
implementation](https://dlmf.nist.gov/10.33). Those functions are accessible
through ``bessel_I_ratio`` and ``bessel_I_ratio_iter``.

# License

This wrapper, including the continued fraction expansion, is licensed under the
[GPLv3](LICENSE.md). The [Fortran algorithm](mod_zbes.f90) published by Masao
Kodama in an ACM journal carries the [ACM license](LICENSE_mod_zbes.md).
