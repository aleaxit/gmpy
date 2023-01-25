# distutils: libraries = gmp mpfr mpc
from gmpy2 cimport *

import_gmpy2()

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

cdef extern from "mpc.h":
    void mpc_init2 (mpc_ptr x, mpfr_prec_t rnd);
    void mpc_clear (mpc_ptr x)

    int mpc_set_si_si (mpc_ptr rop, long int re, long int im, mpc_rnd_t rnd)

def run():
    r"""
    Run all the functions in this module starting with ``test_``
    """
    import traceback, sys
    import test_cython

    failed = 0

    for name in dir(test_cython):
        if not name.startswith('test_'):
            continue
        test = getattr(test_cython, name)
        if not callable(test):
            continue

        sys.stdout.write(name)
        sys.stdout.write('... ')
        sys.stdout.flush()

        try:
            test()
        except Exception:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            sys.stdout.write('FAILED with:\n')
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=2, file=sys.stdout)
            sys.stdout.write('\n')
            failed += 1
        else:
            sys.stdout.write('pass\n')
            sys.stdout.flush()

    if not failed:
        sys.stdout.write('All tests pass!\n')
    else:
        sys.stdout.write('{} failure(s)\n'.format(failed))

def test_mpz():
    cdef mpz x = GMPy_MPZ_New(NULL)
    cdef mpz y = GMPy_MPZ_New(NULL)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpz_set_si(x.z, 3)
    mpz_set_si(y.z, 2)

    assert x != y and y < x and x > y
    assert not ((y == x) or (y >= x))

    # some python operations
    z = x + y + 1
    assert z == z and z == 6 and 6 == z

def test_mpz_cmp():
    cdef mpz z = GMPy_MPZ_New(NULL)

    mpz_set_si(z.z, 3)

    assert z == 3 and 3 == z
    assert z == mpz(3) and mpz(3) == z
    assert z == mpq(3,1) and mpq(3,1) == z
    assert z == 3.0 and 3.0 == z
    assert z == mpfr(3) and mpfr(3) == z
    # Cython issue 1776 prevents testing
    # z == 3 + 0j
    assert z == mpc(3) and mpc(3) == z

def test_mpz_from_mpz():
    cdef mpz_t z
    mpz_init(z)
    mpz_set_si(z, 2)
    cdef mpz a = GMPy_MPZ_From_mpz(z)
    mpz_clear(z)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert a == a and a == 2 and 2 == a

def test_mpq():
    cdef mpq x = GMPy_MPQ_New(NULL)
    cdef mpq y = GMPy_MPQ_New(NULL)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpq_set_si(x.q, 1, 4)
    mpq_set_si(y.q, -3, 7)

    assert x != y

    # some python operations
    z = x + y - 1
    assert z == z and 28 * z == -33 and -33 == z * 28

def test_mpq_cmp():
    cdef mpq q = GMPy_MPQ_New(NULL)

    mpq_set_si(q.q, 1, 1)

    assert q == 1 and 1 == q
    assert q == mpz(1) and mpz(1) == q
    assert q == mpq(1) and mpq(1) == q
    assert q == 1.0 and 1.0 == q
    # Cython issue 1776 prevents testing
    # q == 1 + 0j
    assert q == mpc(1) and mpc(1) == q

def test_mpq_from_mpq():
    cdef mpq_t q
    mpq_init(q)
    mpq_set_si(q, 3, 5)
    cdef mpq a = GMPy_MPQ_From_mpq(q)
    mpq_clear(q)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert a == a and 3 == a * 5 and 5 * a == 3

def test_mpq_from_mpz():
    cdef mpz_t z1, z2
    mpz_init(z1)
    mpz_init(z2)
    mpz_set_si(z1, 51)
    mpz_set_si(z2, 33)
    cdef mpq a = GMPy_MPQ_From_mpz(z1, z2)
    mpz_clear(z1)
    mpz_clear(z2)

    # Check that the refcount is correct
    assert Py_REFCNT(<PyObject *> a) == 1

    assert 33 * a == 51 and 51 == 33 * a

def test_mpfr():
    cdef mpfr x = GMPy_MPFR_New(100, NULL)
    cdef mpfr y = GMPy_MPFR_New(123, NULL)
    cdef mpfr z = GMPy_MPFR_New(0, NULL)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1
    assert Py_REFCNT(<PyObject *> z) == 1

    mpfr_set_si(x.f, 2741, MPFR_RNDN)
    y = x + 1
    z = y + 1

    assert y == y and y == 2742 and 2742 == y
    assert z == z and z == 2743 and 2743 == z

def test_mpfr_cmp():
    cdef mpfr x = GMPy_MPFR_New(53, NULL)

    mpfr_set_si(x.f, 1, MPFR_RNDN)

    assert x == 1 and 1 == x
    assert x == mpz(1) and mpz(1) == x
    assert x == mpq(1) and mpq(1) == x
    assert x == mpfr(1) and mpfr(1) == x
    # Cython issue 1776 prevents testing
    # x == 1 + 0j
    assert x == mpc(1) and mpc(1) == x

def test_mpfr_from_mpfr():
    cdef mpfr_t x
    mpfr_init2(x, 100)
    mpfr_set_si(x, 2341, MPFR_RNDN)
    cdef mpfr a = GMPy_MPFR_From_mpfr(x)
    mpfr_clear(x)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> a) == 1

    assert a == a and a == 2341 and 2341 == a

def test_mpc():
    cdef mpc x = GMPy_MPC_New(56, 56, NULL)
    cdef mpc y = GMPy_MPC_New(56, 56, NULL)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> x) == 1
    assert Py_REFCNT(<PyObject *> y) == 1

    mpc_set_si_si(x.c, 42, 30, MPC_RNDNN)
    y = x * 2

    # work around Cython issue 1776 that prevents testing
    # y == 84 + 60j
    assert y == y and y == mpc(84, 60)

def test_mpc_cmp():
    cdef mpc x = GMPy_MPC_New(0, 0, NULL)
    mpc_set_si_si(x.c, 0, 0, MPC_RNDNN)

    assert x == 0 and 0 == x
    assert x == mpz(0) and mpz(0) == x
    assert x == mpq(0,1) and mpq(0,1) == x
    assert x == 0.0 and 0.0 == x
    # Cython issue 1776 prevents testing
    # x == 0.0j
    assert x == mpc(0,0)

def test_mpc_from_mpc():
    cdef mpc_t x
    mpc_init2(x, 100)
    mpc_set_si_si(x, 2341, 42, MPC_RNDNN)
    cdef mpc a = GMPy_MPC_From_mpc(x)
    mpc_clear(x)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> a) == 1

    # Cython issue 1776 prevents testing
    # a == 2341 + 42j
    assert a == a and a == mpc(2341, 42)

def test_mpc_from_mpfr():
    cdef mpfr_t r
    cdef mpfr_t i

    mpfr_init2(r, 50)
    mpfr_init2(i, 50)

    mpfr_set_si(r, 49, MPFR_RNDN)
    mpfr_set_si(i, 49, MPFR_RNDN)

    cdef mpc c = GMPy_MPC_From_mpfr(r, i)

    mpfr_clear(r)
    mpfr_clear(i)

    # Check that the refcount is appropriate
    assert Py_REFCNT(<PyObject *> c) == 1

    # Cython issue 1776 that prevents testing
    # c == 49 + 49j
    assert c == c and c == mpc(49, 49)

def test_check():
    x = GMPy_MPZ_New(NULL)
    y = GMPy_MPQ_New(NULL)
    z = GMPy_MPFR_New(53, NULL)
    c = GMPy_MPC_New(55, 87, NULL)

    # Checktypes functions with a MPZ_Object
    assert MPZ_Check(x)
    assert not MPQ_Check(x)
    assert not MPC_Check(x)
    assert not MPFR_Check(x)

    # Checktypes functions with a MPQ_Object
    assert not MPZ_Check(y)
    assert MPQ_Check(y)
    assert not MPC_Check(y)
    assert not MPFR_Check(y)

    # Checktypes functions with a MPFR_Object
    assert not MPZ_Check(z)
    assert not MPQ_Check(z)
    assert not MPC_Check(z)
    assert MPFR_Check(z)

    # Checktypes functions with a MPC_Object
    assert not MPZ_Check(c)
    assert not MPQ_Check(c)
    assert MPC_Check(c)
    assert not MPFR_Check(c)

    # Checktypes functions with a python integer
    cdef object i = 5
    assert not MPZ_Check(i)
    assert not MPQ_Check(i)
    assert not MPC_Check(i)
    assert not MPFR_Check(i)

def test_py_check():
    x = mpz()
    y = mpq()
    z = mpfr(53)
    c = mpc(53)

    # Checktypes functions with a MPZ_Object
    assert MPZ_Check(x)
    assert not MPQ_Check(x)
    assert not MPC_Check(x)
    assert not MPFR_Check(x)

    # Checktypes functions with a MPQ_Object
    assert not MPZ_Check(y)
    assert MPQ_Check(y)
    assert not MPC_Check(y)
    assert not MPFR_Check(y)

    # Checktypes functions with a MPFR_Object
    assert not MPZ_Check(z)
    assert not MPQ_Check(z)
    assert not MPC_Check(z)
    assert MPFR_Check(z)

    # Checktypes functions with a MPC_Object
    assert not MPZ_Check(c)
    assert not MPQ_Check(c)
    assert MPC_Check(c)
    assert not MPFR_Check(c)

    # Checktypes functions with a python integer
    cdef object i = 5
    assert not MPZ_Check(i)
    assert not MPQ_Check(i)
    assert not MPC_Check(i)
    assert not MPFR_Check(i)
