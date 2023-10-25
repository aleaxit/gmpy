import math
import numbers
import pickle
import sys
from decimal import Decimal
from fractions import Fraction

import pytest
from hypothesis import example, given, settings
from hypothesis.strategies import fractions, integers
from supportclasses import a, b, c, d, q, z

import gmpy2
from gmpy2 import cmp, cmp_abs, from_binary, mpc, mpfr, mpq, mpz, to_binary


def test_mpz_constructor():
    assert mpq('1/1') == mpq(1,1)
    assert mpq('1_/_1') == mpq(1,1)


def test_mpq_as_integer_ratio():
    assert mpq(2, 3).as_integer_ratio() == (mpz(2), mpz(3))


def test_mpq_numbers_abc():
    assert isinstance(mpq(1, 3), numbers.Rational)


def test_mpq_pickling():
    for proto in range(pickle.HIGHEST_PROTOCOL + 1):
        for x in [mpq(1234, 6789), mpq(-1234, 6789), mpq(0, 1)]:
            assert pickle.loads(pickle.dumps(x, protocol=proto)) == x


def test_mpq_from_float():
    assert mpq.from_float(3.2) == mpq(3602879701896397, 1125899906842624)


def test_mpq_from_Decimal():
    assert mpq(Decimal("5e-3")) == mpq(5, 1000)
    assert mpq(Decimal(1)) == mpq(1)  # issue 327
    assert mpq(Decimal('0.6')) == mpq(3, 5)
    assert mpq.from_decimal(Decimal("5e-3")) == mpq(5, 1000)


def test_mpq_cmp():
    assert cmp(mpq(1,2), mpq(1,2)) == 0
    assert cmp(mpq(-1,2), mpq(1,2)) == -1
    assert cmp(mpq(1,2), mpq(-1,2)) == 1
    assert cmp(0, mpq(1,2)) == -1
    assert cmp(mpq(1,2), 0) == 1

    assert cmp_abs(mpq(1,2), mpq(1,3)) == 1
    assert cmp_abs(mpq(-1,2), mpq(1,3)) == 1
    assert cmp_abs(mpq(1,2), mpq(-1,3)) == 1
    assert cmp_abs(mpq(1,4), mpq(1,3)) == -1
    assert cmp_abs(mpq(-1,4), mpq(1,3)) == -1
    assert cmp_abs(mpq(1,4), mpq(-1,3)) == -1

    assert cmp(mpq(3,2), q) == 0
    assert cmp(q, mpq(3,5)) == 1


def test_mpq_conversion():
    x = mpq(a)
    assert isinstance(x, mpq)
    assert 2*x == 3
    assert mpq(z, z) == mpq(1, 1)

    class Three:
        def __mpz__(self): return mpz(3)

    assert mpq(z, Three()) == mpq(2,3)

    pytest.raises(TypeError, lambda: mpq(b))
    pytest.raises(TypeError, lambda: mpq(c))
    pytest.raises(TypeError, lambda: mpq(d))


def test_mpq_round():
    pytest.raises(TypeError, lambda: mpq(7,3).__round__(4.5))


@settings(max_examples=1000)
@given(integers(), integers(min_value=1))
@example(0, 1)
@example(1, 1)
@example(-1, 2)
@example(123456789123456789, 9876)
def test_mpq_to_from_binary(p, q):
    x = mpq(p,q)
    assert x == from_binary(to_binary(x))


@settings(max_examples=1000)
@given(fractions())
@example(Fraction(0))
@example(Fraction(-1))
@example(Fraction(-2))
@example(Fraction(15432, 125))
@example(Fraction(1, sys.hash_info.modulus))
@example(Fraction(-1, sys.hash_info.modulus))
def test_mpq_hash(q):
    assert hash(mpq(q)) == hash(q)


def test_mpq_digits():
    q = mpq(2/3)

    assert q.digits() == '6004799503160661/9007199254740992'
    assert q.digits(16) == '0x15555555555555/0x20000000000000'

    pytest.raises(TypeError, lambda: q.digits(16, 5))
    pytest.raises(ValueError, lambda: q.digits(0))


def test_mpq_sub():
    assert mpq(1,2) - Fraction(3,2) == mpq(-1,1)
    assert Fraction(1,2) - mpq(3,2) == mpq(-1,1)
    assert mpq(1,2) - mpq(3,2) == mpq(-1,1)
    assert mpq(1,2) - 0 == mpq(1,2)
    assert mpq(1,2) - mpz(1) == mpq(-1,2)
    assert mpq(1,2) + (-1) == mpq(-1,2)
    assert 1 - mpq(1,2) == mpq(1,2)
    assert mpz(1) - mpq(1,2) == mpq(1,2)
    assert mpq(1,2) - mpz(1) == mpq(-1,2)
    assert mpq(1,1) - mpc(0) == mpc('1.0+0.0j')
    assert mpc(0) - mpq(1,1) == mpc('-1.0+0.0j')
    assert mpq(1,2) - z == mpq(-3,2)
    assert mpq(1,2) - q ==mpq(-1,1)

    ctx = gmpy2.context()

    assert ctx.sub(mpq(1,2), mpq(3,2)) == mpq(-1,1)
    assert ctx.sub(mpq(1,2), Fraction(3,2)) == mpq(-1,1)
    assert ctx.sub(Fraction(1,2), mpq(3,2)) == mpq(-1,1)
    assert ctx.sub(Fraction(1,2), Fraction(3,2)) == mpq(-1,1)

    pytest.raises(TypeError, lambda: ctx.sub(1,'a'))
    pytest.raises(TypeError, lambda: mpq(1,2) - 'a')
    pytest.raises(TypeError, lambda: 'a' - mpq(1,2))


def test_mpq_mul():
    assert mpq(1,2) * Fraction(3,2) == mpq(3,4)
    assert Fraction(1,2) * mpq(3,2) == mpq(3,4)
    assert mpq(1,2) * mpq(3,2) == mpq(3,4)
    assert mpq(1,2) * 0 == mpq(0,1)
    assert mpq(1,2) * mpz(1) == mpq(1,2)
    assert mpq(1,2) * (-1) == mpq(-1,2)
    assert mpq(1,1) * mpc(1,0) == mpc('1.0+0.0j')
    assert mpc(1,0) * mpq(1,1) == mpc('1.0+0.0j')
    assert mpq(1,2) * z == mpq(1,1)
    assert mpq(1,2) * q == mpq(3,4)

    ctx = gmpy2.context()

    assert ctx.mul(mpq(1,2), mpq(3,2)) == mpq(3,4)
    assert ctx.mul(mpq(1,2), Fraction(3,2)) == mpq(3,4)
    assert ctx.mul(Fraction(1,2), mpq(3,2)) == mpq(3,4)
    assert ctx.mul(Fraction(1,2), Fraction(3,2)) == mpq(3,4)

    pytest.raises(TypeError, lambda: ctx.mul(1,'a'))
    pytest.raises(TypeError, lambda: mpq(1,2) * 'a')
    pytest.raises(TypeError, lambda: 'a' * mpq(1,2))


def test_mpq_divmod():
    pytest.raises(TypeError, lambda: divmod(mpq(1,2),'a'))

    ctx = gmpy2.ieee(64)
    gmpy2.set_context(ctx)

    assert ctx.divmod(mpq(3,2),mpq(3,7)) == (mpz(3), mpq(3,14))

    pytest.raises(TypeError, lambda: divmod(mpq(1,2), mpc(1,2)))


def test_mpq_attributes():
    q = mpq('4/5')
    pyq = Fraction(4, 5)

    assert q.numerator == mpz(4)
    assert q.denominator == mpz(5)
    assert gmpy2.numer(q) == mpz(4)
    assert gmpy2.denom(q) == mpz(5)

    pytest.raises(TypeError, lambda: gmpy2.numer(6.2))
    pytest.raises(TypeError, lambda: gmpy2.denom(5.6))
    pytest.raises(TypeError, lambda: gmpy2.denom(mpfr(5)))

    assert gmpy2.numer(pyq) == mpz(4)
    assert gmpy2.denom(pyq) == mpz(5)


def test_mpq_floor():
    assert math.floor(mpq('7/2')) == mpz(3)


def test_mpq_ceil():
    assert math.ceil(mpq('7/2')) == mpz(4)


def test_mpq_trunc():
    assert math.trunc(mpq('7/2')) == mpz(3)


def test_mpq_round():
    q = mpq('4/5')

    assert round(mpq('7/2')) == mpz(4)
    assert round(q, 4) == mpq(4,5)

    pytest.raises(TypeError, lambda: round(q, 4, 2))


def test_mpq_not():
    q = mpq('4/5')

    assert q
    assert not mpq('0/5')


def test_mpq_qdiv():
    q = mpq('4/5')
    pyq = Fraction(4, 5)

    assert gmpy2.qdiv(q) == mpq(4,5)
    assert gmpy2.qdiv(pyq) == mpq(4,5)
    assert gmpy2.qdiv(5) == mpz(5)

    pytest.raises(TypeError, lambda: gmpy2.qdiv(mpc(4, 5)))
    pytest.raises(TypeError, lambda: gmpy2.qdiv(4, 5, 4))
    pytest.raises(TypeError, lambda: gmpy2.qdiv(4, 5.6))

    assert gmpy2.qdiv(q, 2) == mpq(2,5)
    assert gmpy2.qdiv(10, q) == mpq(25,2)
    assert gmpy2.qdiv(1) == mpz(1)


def test_issue_334():
    x = mpq(3,2)
    y = mpq(x,2)

    assert x == mpq(3,2)
    assert y == mpq(3,4)
    assert id(x) is not id(y)
