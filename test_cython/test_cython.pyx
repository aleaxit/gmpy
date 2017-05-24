# distutils: libraries = gmp
from __future__ import print_function
from gmpy2 cimport *

import_gmpy2()

import sys

cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass

    long Py_REFCNT(PyObject *)

cdef extern from "gmp.h":
    void mpz_set_si(mpz_t, long)
    void mpq_set_si(mpq_t, long, long)

def test_mpz():
    sys.stdout.write('test mpz... ')
    sys.stdout.flush()

    x = MPZ_New()
    y = MPZ_New()

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpz_set_si(MPZ(<MPZ_Object *> x), 3)
    mpz_set_si(MPZ(<MPZ_Object *> y), 2)

    # some python operations
    z = x + y + 1
    assert z == 6

    sys.stdout.write('done\n')
    sys.stdout.flush()


def test_mpq():
    sys.stdout.write('test mpq... ')
    sys.stdout.flush()

    x = MPQ_New()
    y = MPQ_New()

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpq_set_si(MPQ(<MPQ_Object *> x), 1, 4)
    mpq_set_si(MPQ(<MPQ_Object *> y), -3, 7)

    # some python operations
    z = x + y - 1
    assert 28 * z == -33

    sys.stdout.write('done\n')
    sys.stdout.flush()


def test_check():
    sys.stdout.write('test check... ')
    sys.stdout.flush()

    x = MPZ_New()
    y = MPQ_New()

    # Checktypes functions with a MPZ_Object
    assert MPZ_Check(<PyObject *>  x)
    assert not MPQ_Check(<PyObject *> x)
    assert not MPC_Check(<PyObject *> x)
    assert not MPFR_Check(<PyObject *> x)

    # Checktypes functions with a MPQ_Object
    assert not MPZ_Check(<PyObject *> y)
    assert MPQ_Check(<PyObject *> y)
    assert not MPC_Check(<PyObject *> y)
    assert not MPFR_Check(<PyObject *> y)

    # Checktypes functions with a python integer
    cdef object i = 5
    assert not MPZ_Check(<PyObject *> i)
    assert not MPQ_Check(<PyObject *> i)
    assert not MPC_Check(<PyObject *> i)
    assert not MPFR_Check(<PyObject *> i)

    sys.stdout.write('done\n')
    sys.stdout.flush()

test_mpz()
test_mpq()
test_check()
