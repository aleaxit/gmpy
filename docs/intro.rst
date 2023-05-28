Introduction
============

gmpy2 is a C-coded Python extension module that supports multiple-precision
arithmetic.  It is the successor to the original gmpy module (supported only
the GMP library). gmpy2 adds support for the MPFR (correctly rounded real
floating-point arithmetic) and MPC (correctly rounded complex floating-point
arithmetic) libraries.

The following libraries are supported:

* GMP for integer and rational arithmetic (https://gmplib.org).
* MPFR for correctly rounded real floating-point arithmetic
  (https://www.mpfr.org).
* MPC for correctly rounded complex floating-point arithmetic
  (https://mpc.multiprecision.org).
* Generalized Lucas sequences and primality tests are based on the following
  code:

      - mpz_lucas: https://sourceforge.net/projects/mpzlucas/
      - mpz_prp: https://sourceforge.net/projects/mpzprp/

gmpy2 Versions
--------------

gmpy2 2.1.x are --- last releases that will support Python versions <= 3.7.

Installation
============

Pre-compiled binary wheels are available on PyPI for Linux, MacOS, and Windows,
you can install latest release of the gmpy2 with pip::

    pip install gmpy2

or some specific version with::

    pip install gmpy2==2.1.5
