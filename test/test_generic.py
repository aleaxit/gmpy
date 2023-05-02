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
