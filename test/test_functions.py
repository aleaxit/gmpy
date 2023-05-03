import pytest

from gmpy2 import (root, rootn, zero, mpz, mpq, mpfr, mpc, is_nan, maxnum,
                   minnum)


def test_root():
    assert root(zero(1), 2) == mpfr('0.0')
    assert root(zero(-1), 2) == mpfr('-0.0')
    assert root(zero(-1), 3) == mpfr('-0.0')
    assert root(zero(-1), 4) == mpfr('-0.0')
    assert root(zero(-1), 5) == mpfr('-0.0')
    assert root(zero(-1), 6) == mpfr('-0.0')
    assert root(2, 2) == mpfr('1.4142135623730951')
    assert root(mpz(2), 2) == mpfr('1.4142135623730951')
    assert root(mpq(2), 2) == mpfr('1.4142135623730951')
    assert root(mpfr(2), 2) == mpfr('1.4142135623730951')
    pytest.raises(TypeError, lambda: root(mpc(2), 2))
    assert is_nan(root(-2, 2))
    pytest.raises(OverflowError, lambda: root(2, -2))
    pytest.raises(TypeError, lambda: root(2, 0.5))


def test_rootn():
    assert rootn(zero(1), 2) == mpfr('0.0')
    assert rootn(zero(-1), 2) == mpfr('0.0')
    assert rootn(zero(-1), 3) == mpfr('-0.0')
    assert rootn(zero(-1), 4) == mpfr('0.0')
    assert rootn(zero(-1), 5) == mpfr('-0.0')
    assert rootn(zero(-1), 6) == mpfr('0.0')


def test_maxnum():
    a = mpfr("12.34")
    b = mpfr("45.67")
    nan = mpfr("nan")
    inf = mpfr("inf")

    assert maxnum(a, b) == mpfr('45.670000000000002')
    assert maxnum(b, a) == mpfr('45.670000000000002')
    assert maxnum(a, -b) == mpfr('12.34')
    assert maxnum(a, 123456) == mpfr('123456.0')
    assert maxnum(12345678901234567890, a) == mpfr('1.2345678901234567e+19')
    assert maxnum(0, -1) == mpfr('0.0')
    assert maxnum(1, inf) == mpfr('inf')
    assert maxnum(1, -inf) == mpfr('1.0')
    assert maxnum(nan, a) == mpfr('12.34')
    assert maxnum(a, nan) == mpfr('12.34')
    assert maxnum(nan, inf) == mpfr('inf')
    assert maxnum(nan, -inf) == mpfr('-inf')
    assert is_nan(maxnum(nan, nan))


def test_minnum():
    a = mpfr("12.34")
    b = mpfr("45.67")
    nan = mpfr("nan")
    inf = mpfr("inf")
    minf = mpfr("-inf")

    assert minnum(a, b) == mpfr('12.34')
    assert minnum(b, a) == mpfr('12.34')
    assert minnum(1, inf) == mpfr('1.0')
    assert minnum(minf, a) == mpfr('-inf')
    assert minnum(nan, inf) == mpfr('inf')
    assert is_nan(minnum(nan, nan))
