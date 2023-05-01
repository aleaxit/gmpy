import pytest

from gmpy2 import root, rootn, zero, mpz, mpq, mpfr, mpc, is_nan


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
