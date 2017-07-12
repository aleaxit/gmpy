Cython usage
============

The gmpy2 module provides a C-API that can be conveniently used from Cython.
All types and functions are declared in the header gmpy2.pxd that is installed
automatically in your Python path together with the library.

Initialization
--------------

In order to use the C-API you need to make one call to the function **void import_gmpy2(void)**.

Types
-----

The mpz, mpq, mpfr and mpc Python types correspond to the C struct

**struct MPZ_Object ***

**struct MPQ_Object ***

**struct MPFR_Object ***

**struct MPC_Object ***

Fast type checking can be done with the following C functions

**bint MPZ_Check(PyObject * obj)**
    equivalent to **isinstance(obj, mpz)**

**bint MPQ_Check(PyObject * obj)**
    equivalent to **isinstance(obj, mpq)**

**bint MPFR_Check(PyObject * obj)**
    equivalent to **isinstance(obj, mpfr)**

**bint MPC_Check(PyObject * obj)**
    equivalent to **isinstance(obj, mpc)**

Object creation
---------------

To create a new gmpy2 types there are four basic functions

**object MPZ_New()**
    create a new mpz object

**object MPQ_New()**
    create a new mpq object

**object MPFR_New(mpfr_prec_t prec)**
    create a new mpfr object with given precision prec

**object  MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec)**
    create a new mpc object with given precisions rprec and iprec of
    respectively real and imaginary parts

The gmpy2.pxd header also provides convenience macro to wrap a (copy of) a mpz_t, mpq_t, mpfr_t
or a mpc_t object into the corresponding gmpy2 type.

**object GMPy_MPZ_From_mpz(mpz_srcptr z)**
    return a new mpz object with a given mpz_t value z

**object GMPy_MPQ_From_mpq(mpq_srcptr q)**
    return a new mpq object from a given mpq_t value q

**object GMPy_MPQ_From_mpz(mpz_srcptr num, mpz_srcptr den)**
    return a new mpq object with a given mpz_t numerator num and mpz_t denominator den

**object GMPy_MPFR_From_mpfr(mpfr_srcptr x)**
    return a new mpfr object with a given mpfr_t value x

**object GMPy_MPC_From_mpc(mpc_srcptr c)**
    return a new mpc object with a given mpc_t value c

**object GMPy_MPC_From_mpfr(mpfr_srcptr re, mpfr_srcptr im)**
    return a new mpc object with a given mpfr_t real part re and mpfr_t imaginary part im

Access to the underlying C type
--------------------------------

Each of the gmpy2 objects has a field corresponding to the underlying C
type. The following functions give access to this field

**mpz_t MPZ(MPZ_Object *)**

**mpq_t MPQ(MPQ_OBJECT *)**

**mpfr_t MPFR(MPFR_Object *)**

**mpc_t MPC(MPC_Object *)**

Compilation
------------

The header gmpy2.pxd as well as the C header gmpy2.h from which it depends
are installed in the Python path. In order to make Cython and the C compiler aware
of the existence of these files, the Python path should be part of the include
directories.

Recall that **import_gmpy2()** needs to be called *before* any other function of
the C-API.

Here is a minimal example of a Cython file test_gmpy2.pyx

::

    "A minimal setup.py for compiling test_gmpy2.pyx"

    from distutils.core import setup
    from distutils.extension import Extension
    from Cython.Build import cythonize
    import sys

    ext = Extension("test_gmpy2", ["test_gmpy2.pyx"], include_dirs=sys.path, libraries=['gmp', 'mpfr', 'mpc'])

    setup(
        name="cython_gmpy_test",
        ext_modules=cythonize([ext], include_path=sys.path)
    )

The corresponding setup.py is given below.

::

    "A minimal cython file test_gmpy2.pyx"

    from gmpy2 cimport *

    cdef extern from "gmp.h":
        void mpz_set_si(mpz_t, long)

    import_gmpy2()   # needed to initialize the C-API

    z = MPZ_New()
    mpz_set_si(MPZ(<MPZ_Object *> z), -7)

    print(z + 3)


With these two files in the same repository, you should be able to compile your
module using

::

    $ python setup.py build_ext --inplace

For more about compilation and installation of cython files and extension
modules, please refer to the official documentation of Cython and distutils.
