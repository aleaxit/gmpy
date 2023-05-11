from fractions import Fraction

import pytest

from gmpy2 import mpz, mpq, mpfr, mpc, get_context, square, add
from supportclasses import z, q, r, cx


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


def test_square():
    z = mpz(2)
    assert square(z) == mpz(4)
    assert square(z) == z * z

    q = mpq(2,3)
    assert square(q) == mpq(4,9)
    assert square(q) == q * q

    r = mpfr(5.3)
    assert square(r) == mpfr('28.09')
    assert square(r) == r * r

    c = mpc(2,3)
    assert square(c) == mpc('-5.0+12.0j')
    assert square(c) == c * c
    assert square(2) == square(z)
    assert square(Fraction(2,3)) == square(q)
    assert square(5.3) == square(r)
    assert square(complex(2,3)) == square(c)

    pytest.raises(TypeError, lambda: square('invalid'))


def test_add():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890

    assert a+1 == mpz(124)
    assert a+(-1) == mpz(122)
    assert 1+a == mpz(124)
    assert (-1)+a == mpz(122)

    assert a+c == mpz(12345678901234568013)
    assert c+a == mpz(12345678901234568013)
    assert a+(-c) == mpz(-12345678901234567767)
    assert (-c)+a == mpz(-12345678901234567767)

    assert a+b == mpz(579)
    assert b+a == mpz(579)
    assert a+(-b) == mpz(-333)
    assert (-b)+a == mpz(-333)
    assert a+z == mpz(125)

    assert add(a,b) == a+b
    assert add(c,c) == c+c
    assert add(1, 1) == mpz(2)
    assert add(mpfr(1), mpfr(1)) == mpfr('2.0')
    assert add(a, 1) == mpz(124)
    assert add(1, a) == mpz(124)
    assert add(a, mpq(0)) == mpq(123,1)
    assert add(a, mpfr(0)) == mpfr('123.0')
    assert add(a, mpc(0)) == mpc('123.0+0.0j')

    pytest.raises(TypeError, lambda: add(1))
    pytest.raises(TypeError, lambda: add(1,2,3))

    assert add(1,1) == mpz(2)

    pytest.raises(TypeError, lambda: a+'b')
    pytest.raises(TypeError, lambda: 'b'+a)

    assert mpq(1,2) + Fraction(3,2) == mpq(2,1)
    assert Fraction(1,2) + mpq(3,2) == mpq(2,1)
    assert mpq(1,2) + mpq(3,2) == mpq(2,1)
    assert mpq(1,2) + 0 == mpq(1,2)
    assert mpq(1,2) + mpz(1) == mpq(3,2)
    assert mpq(1,2) + (-1) == mpq(-1,2)
    assert 1 + mpq(1,2) == mpq(3,2)
    assert mpz(1) + mpq(1,2) == mpq(3,2)
    assert mpq(1,2) + mpz(1) == mpq(3,2)
    assert mpq(1,1) + mpc(0) == mpc('1.0+0.0j')
    assert mpc(0) + mpq(1,1) == mpc('1.0+0.0j')
    assert mpq(1,2) + z == mpq(5,2)
    assert mpq(1,2) + q == mpq(2,1)

    assert add(mpq(1,2), mpq(3,2)) == mpq(2,1)
    assert add(mpq(1,2), Fraction(3,2)) == mpq(2,1)
    assert add(Fraction(1,2), mpq(3,2)) == mpq(2,1)
    assert add(Fraction(1,2), Fraction(3,2)) == mpq(2,1)

    pytest.raises(TypeError, lambda: add(1, 'a'))
    pytest.raises(TypeError, lambda: mpq(1,2) + 'a')
    pytest.raises(TypeError, lambda: 'a' + mpq(1,2))

    assert mpfr(1) + 1 == mpfr('2.0')
    assert 1 + mpfr(1) == mpfr('2.0')
    assert mpfr(1) + mpz(1) == mpfr('2.0')
    assert mpz(1) + mpfr(1) == mpfr('2.0')
    assert mpfr(1) + mpfr(1) == mpfr('2.0')
    assert mpfr(1) + mpq(1,1) == mpfr('2.0')
    assert mpq(1,1) + mpfr(1) == mpfr('2.0')
    assert mpfr(1) + Fraction(1,1) == mpfr('2.0')
    assert Fraction(1,1) + mpfr(1) == mpfr('2.0')
    assert mpfr(1) + 1.0 == mpfr('2.0')
    assert 1.0 + mpfr(1) == mpfr('2.0')
    assert mpfr(0) + (1 << 100) == mpfr('1p100', base=2)
    assert (1 << 100) + mpfr(0) == mpfr('1p100', base=2)
    assert mpfr(1) + z == mpfr('3.0')
    assert mpfr(0.5) + q == mpfr('2.0')
    assert mpfr(1.5) + r == mpfr('3.0')

    pytest.raises(TypeError, lambda: mpc(1,2) + 'a')

    assert mpfr(1) + mpc(1,2) == mpc('2.0+2.0j')
    assert mpc(1,2) + mpfr(1) == mpc('2.0+2.0j')
    assert mpc(1,2) + 1+0j == mpc('2.0+2.0j')
    assert 1+0j + mpc(1,2) == mpc('2.0+2.0j')
    assert mpc(1,2) + cx == mpc('43.0+69.0j')
    assert mpc(1,2) + r == mpc('2.5+2.0j')
    assert mpc(1,2) + q == mpc('2.5+2.0j')
    assert mpc(1,2) + z == mpc('3.0+2.0j')
