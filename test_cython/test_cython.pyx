# distutils: libraries = gmp
from __future__ import print_function
from gmpy2 cimport *

import_gmpy2()

cdef extern from "gmp.h":
    void mpz_set_si(mpz_t, long)
    void mpq_set_si(mpq_t, long, long)
    void mpz_add(mpz_t, mpz_t, mpz_t)

cdef x = <object> GMPy_MPZ_New(NULL)
cdef y = <object> GMPy_MPZ_New(NULL)

mpz_set_si(MPZ(<MPZ_Object *> x), 3)
mpz_set_si(MPZ(<MPZ_Object *> y), 2)

z = x + y + 1  # python operation!
assert z == 6

# Checktypes functions with a MPZ_Object
assert MPZ_Check(<PyObject *>x)
assert not MPQ_Check(<PyObject *>x)
assert not MPC_Check(<PyObject *>x)
assert not MPFR_Check(<PyObject *>x)

# Checktypes functions with a python integer
i=5
assert not MPZ_Check(<PyObject *>i)
assert not MPQ_Check(<PyObject *>i)
assert not MPC_Check(<PyObject *>i)
assert not MPFR_Check(<PyObject *>i)

x = <object> GMPy_MPQ_New(NULL)
y = <object> GMPy_MPQ_New(NULL)

mpq_set_si(MPQ(<MPQ_Object *> x), 1, 4)
mpq_set_si(MPQ(<MPQ_Object *> y), -3, 7)

z = x + y - 1 # python operation!
assert 28 * z == -33

# Checktypes functions with a MPZ_Object
assert not MPZ_Check(<PyObject *>x)
assert MPQ_Check(<PyObject *>x)
assert not MPC_Check(<PyObject *>x)
assert not MPFR_Check(<PyObject *>x)
