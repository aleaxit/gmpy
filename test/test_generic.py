from fractions import Fraction

import pytest

from gmpy2 import mpz, mpq, mpfr, mpc, get_context


def test_minus():
    ctx = get_context()
    assert ctx.minus(mpz(5)) == mpz(-5)
    assert ctx.minus(mpq(4,5)) == mpq(-4,5)
    assert ctx.minus(mpfr(4,5)) == mpfr('-4.0')
    assert ctx.minus(mpfr('inf')) == mpfr('-inf')
    assert ctx.minus(mpc(15,3)) == mpc('-15.0-3.0j')
    assert ctx.minus(65) == mpz(-65)
    assert ctx.minus(5.5) == mpfr('-5.5')
    assert ctx.minus(Fraction(2,3)) == mpq(-2,3)
    assert ctx.minus(complex(15,3)) == mpc('-15.0-3.0j')
    assert -mpc(5,5) == mpc('-5.0-5.0j')

    pytest.raises(TypeError, lambda: ctx.minus('invalid'))
    pytest.raises(TypeError, lambda: ctx.minus())


def test_plus():
    ctx = get_context()
    assert ctx.plus(5) == mpz(5)
    assert ctx.plus(-5) == mpz(-5)
    assert ctx.plus(Fraction(4,5)) == mpq(4,5)
    assert ctx.plus(4.5) == mpfr('4.5')
    assert ctx.plus(complex(5.2,5)) == mpc('5.2000000000000002+5.0j')
    assert ctx.plus(mpz(421)) == mpz(421)

    pytest.raises(TypeError, lambda: ctx.plus('invalid'))
    pytest.raises(TypeError, lambda: ctx.plus())

    assert + mpz(421) == mpz(421)
    assert + mpq('4/5') == mpq(4,5)
    assert + mpfr('inf') == mpfr('inf')
    assert + mpc(65.0, 45) == mpc('65.0+45.0j')
