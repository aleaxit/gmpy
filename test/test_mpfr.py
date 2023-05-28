from decimal import Decimal

import pytest
from hypothesis import given, example, settings
from hypothesis.strategies import floats
import mpmath

import gmpy2
from gmpy2 import (gamma_inc, mpfr, cmp, cmp_abs, zero, nan, mpz, mpq,
                   to_binary, from_binary, is_nan, random_state,
                   mpfr_grandom, mpfr_nrandom)
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


def test_mpfr_hash():
    assert hash(mpfr('123.456')) == hash(float('123.456'))
    assert hash(mpfr('123.5')) == hash(float('123.5'))
    assert hash(mpfr('0')) == hash(float('0'))
    assert hash(mpfr('1')) == hash(float('1'))
    assert hash(mpfr('2')) == hash(float('2'))
    assert hash(mpfr('-1')) == hash(float('-1'))
    assert hash(mpfr('Inf')) == hash(float('Inf'))
    assert hash(mpfr('-Inf')) == hash(float('-Inf'))
    assert hash(mpfr('-0')) == hash(float('-0'))
    assert hash(mpfr('123.456')) != hash(Decimal('123.456'))
    assert hash(mpfr('123.5')) == hash(Decimal('123.5'))


@given(floats())
@example(0.0)
@example(1.0)
@example(-1.0)
@example(+float('inf'))
@example(-float('inf'))
@example(float('nan'))
@example(1.345)
def test_mpfr_to_from_binary_bulk(r):
    x = mpfr(r)
    y = from_binary(to_binary(x))
    assert x == y or (is_nan(x) and is_nan(y))


def test_mpfr_to_from_binary():
    x = mpfr("1.345e1000")
    assert x==from_binary(to_binary(x))
    x = gmpy2.const_pi()
    assert x.rc == -1
    y = from_binary(to_binary(x))
    assert x == y and y.rc == -1
    -1
    with gmpy2.local_context() as ctx:
        ctx.precision = 100
        x = gmpy2.const_pi()
        assert x == from_binary(to_binary(x))
        ctx.precision = 200
        x = mpfr(gmpy2.const_pi())
        assert x == from_binary(to_binary(x))
        x = gmpy2.const_pi()
        ctx.precision = 300
        x = from_binary(to_binary(x))
        assert x.precision == 200


def test_mpfr_random():
    assert gmpy2.mpfr_random(random_state(42)) == mpfr('0.93002690534702315')


def test_mpfr_grandom():
    assert mpfr_grandom(random_state(42)) == (mpfr('-0.32898912492644183'),
                                              mpfr('0.03656576719642516'))


def test_mpfr_nrandom():
    assert mpfr_nrandom(random_state(42)) == mpfr('-0.32898912492644183')


def test_mpfr_mpmath():
    a, b, c, d = '1.1', '-1.1', '-3.14', '0'
    assert mpfr(a)._mpf_ == (0, mpz(4953959590107546), -52, 53)
    assert mpmath.mpf(mpfr(a)) == mpmath.mpf(a)
    assert mpfr(b)._mpf_ == (1, mpz(4953959590107546), -52, 53)
    assert mpmath.mpf(mpfr(b)) == mpmath.mpf(b)
    assert mpfr(c, precision=10)._mpf_ == (1, mpz(804), -8, 10)
    assert mpmath.mpf(mpfr(c, precision=10), prec=10) == mpmath.mpf(c, prec=10)
    assert mpfr(d)._mpf_ == (0, mpz(0), 1, 1)
