# distutils: libraries = gmp mpfr
from __future__ import print_function
from gmpy2 cimport *

import_gmpy2()

import sys

cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass

    long Py_REFCNT(PyObject *)

cdef extern from "gmp.h":
    void mpz_init (mpz_t integer)
    void mpz_clear (mpz_t integer)

    void mpq_init (mpq_t dest_rational)
    void mpq_clear (mpq_t rational_number)

    void mpz_set_si(mpz_t, long)
    void mpq_set_si(mpq_t, long, long)

cdef extern from "mpfr.h":
    void mpfr_init2 (mpfr_t x, mpfr_prec_t prec)
    void mpfr_clear (mpfr_t x)

    int mpfr_set_si (mpfr_t rop, long int op, mpfr_rnd_t rnd)

def test_mpz():
    sys.stdout.write('test mpz... ')
    sys.stdout.flush()

    x = MPZ_New()
    y = MPZ_New()

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpz_set_si(MPZ(<MPZ_Object *> x), 3)
    mpz_set_si(MPZ(<MPZ_Object *> y), 2)

    # some python operations
    z = x + y + 1
    assert z == 6

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpz_from_mpz():
    sys.stdout.write('test mpz from mpz...')
    sys.stdout.flush()

    cdef mpz_t z
    mpz_init(z)
    mpz_set_si(z, 2)
    a = GMPy_MPZ_From_mpz(z)
    mpz_clear(z)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert a == 2

    del a

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpq():
    sys.stdout.write('test mpq... ')
    sys.stdout.flush()

    x = MPQ_New()
    y = MPQ_New()

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpq_set_si(MPQ(<MPQ_Object *> x), 1, 4)
    mpq_set_si(MPQ(<MPQ_Object *> y), -3, 7)

    # some python operations
    z = x + y - 1
    assert 28 * z == -33

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpq_from_mpq():
    sys.stdout.write('test mpq from mpq...')
    sys.stdout.flush()

    cdef mpq_t q
    mpq_init(q)
    mpq_set_si(q, 3, 5)
    a = GMPy_MPQ_From_mpq(q)
    mpq_clear(q)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert 5 * a == 3

    del a

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpq_from_mpz():
    sys.stdout.write('test mpq from mpz...')
    sys.stdout.flush()

    cdef mpz_t z1, z2
    mpz_init(z1)
    mpz_init(z2)
    mpz_set_si(z1, 51)
    mpz_set_si(z2, 33)
    a = GMPy_MPQ_From_mpz(z1, z2)
    mpz_clear(z1)
    mpz_clear(z2)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert 33 * a == 51

    del a

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpfr():
    sys.stdout.write('test mpfr... ')
    sys.stdout.flush()

    x = MPFR_New(100)
    y = MPFR_New(123)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpfr_set_si(MPFR(<MPFR_Object *> x), 2741, MPFR_RNDN)
    y = x + 1
    assert y == 2742

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_mpfr_from_mpfr():
    sys.stdout.write('test mpfr from mpfr... ')
    sys.stdout.flush()

    cdef mpfr_t x
    mpfr_init2(x, 100)
    mpfr_set_si(x, 2341, MPFR_RNDN)
    a = GMPy_MPFR_From_mpfr(x)
    mpfr_clear(x)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> a) == 1

    assert a == 2341

    del a

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_check():
    sys.stdout.write('test check... ')
    sys.stdout.flush()

    x = MPZ_New()
    y = MPQ_New()
    z = MPFR_New(53)

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

    # Checktypes functions with a MPFR_Object
    assert not MPZ_Check(<PyObject *>  z)
    assert not MPQ_Check(<PyObject *> z)
    assert not MPC_Check(<PyObject *> z)
    assert MPFR_Check(<PyObject *> z)

    # Checktypes functions with a python integer
    cdef object i = 5
    assert not MPZ_Check(<PyObject *> i)
    assert not MPQ_Check(<PyObject *> i)
    assert not MPC_Check(<PyObject *> i)
    assert not MPFR_Check(<PyObject *> i)

    sys.stdout.write('done\n')
    sys.stdout.flush()

def test_py_check():
    sys.stdout.write('test py check... ')
    sys.stdout.flush()

    from gmpy2 import mpz, mpq, mpfr
    x = mpz()
    y = mpq()
    z = mpfr(53)

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

    # Checktypes functions with a MPFR_Object
    assert not MPZ_Check(<PyObject *>  z)
    assert not MPQ_Check(<PyObject *> z)
    assert not MPC_Check(<PyObject *> z)
    assert MPFR_Check(<PyObject *> z)

    # Checktypes functions with a python integer
    cdef object i = 5
    assert not MPZ_Check(<PyObject *> i)
    assert not MPQ_Check(<PyObject *> i)
    assert not MPC_Check(<PyObject *> i)
    assert not MPFR_Check(<PyObject *> i)

    sys.stdout.write('done\n')
    sys.stdout.flush()

test_mpz()
test_mpz_from_mpz()
test_mpq_from_mpz()
test_mpq()
test_mpq_from_mpq()
test_mpfr()
test_mpfr_from_mpfr()
test_check()
test_py_check()
