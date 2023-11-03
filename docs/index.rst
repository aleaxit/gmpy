Welcome to gmpy2's documentation!
=================================

gmpy2 is a C-coded Python extension module that supports multiple-precision
arithmetic.  It is the successor to the original gmpy module (supported only
the GMP library). gmpy2 adds support for the MPFR (correctly rounded real
floating-point arithmetic) and MPC (correctly rounded complex floating-point
arithmetic) libraries.

The following libraries are supported:

* GMP for integer and rational arithmetic (https://gmplib.org)
* MPFR (https://www.mpfr.org)
* MPC (https://www.multiprecision.org/mpc/)
* Generalized Lucas sequences and primality tests are based on the following
  code:

      * mpz_lucas: https://sourceforge.net/projects/mpzlucas/
      * mpz_prp: https://sourceforge.net/projects/mpzprp/

Contents
--------

.. toctree::
   :maxdepth: 2

   overview
   install
   tutorial
   mpz
   advmpz
   mpq
   contexts
   exceptions
   mpfr
   mpc
   generic
   misc
   cython
   conversion
   history

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`
