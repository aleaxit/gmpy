import pytest

import gmpy2
from gmpy2 import gamma_inc, mpfr, cmp, cmp_abs, zero, nan, mpz, mpq
from supportclasses import a, b, c, d, q, r


def test_mpfr_gamma_inc():
    assert gamma_inc(1, 1) == mpfr('0.36787944117144233')
    assert gamma_inc(1, 0) == mpfr('1.0')
    assert gamma_inc(0, 1) == mpfr('0.21938393439552029')


def test_mpfr_cmp():
    assert cmp(mpfr(0), mpfr(0)) == 0
    assert cmp(mpfr(0), mpz(0)) == 0
    assert cmp(mpfr(0), mpq(0,1)) == 0
    assert cmp(zero(-1), zero(-1)) == 0
    assert cmp(zero(1), zero(-1)) == 0
    assert cmp(zero(-1), zero(1)) == 0
    assert cmp(mpfr(1), mpfr(0)) == 1
    assert cmp(mpfr(1), mpz(0)) == 1
    assert cmp(mpfr(1), mpq(0,1)) == 1
    assert cmp(mpfr(-1), mpfr(0)) == -1
    assert cmp(mpfr(-1), mpz(0)) == -1
    assert cmp(mpfr(-1), mpq(0,1)) == -1
    assert cmp(nan(), mpfr(0)) == 0
    assert cmp(nan(), mpz(0)) == 0
    assert cmp(nan(), mpq(0,1)) == 0

    gmpy2.get_context().clear_flags()

    assert cmp(nan(), 1) == 0
    assert gmpy2.get_context().erange is True

    assert cmp_abs(mpfr(-1), mpfr(0)) == 1
    assert cmp_abs(mpfr(-1), mpz(0)) == 1
    assert cmp_abs(mpfr(-1), mpq(0,1)) == 1
    assert cmp_abs(mpfr(0), mpfr(-1)) == -1
    assert cmp_abs(mpz(0), mpfr(-1)) == -1
    assert cmp_abs(mpq(0,1), mpfr(-1)) == -1

    assert cmp(mpfr(1.5), q) == 0
    assert cmp(r, mpfr(1.5)) == 0


def test_mpfr_conversion():
    x = mpfr(a)
    assert isinstance(x, mpfr)
    assert x == 1.5
    pytest.raises(TypeError, lambda: mpfr(b))
    pytest.raises(TypeError, lambda: mpfr(c))
    pytest.raises(TypeError, lambda: mpfr(d))
