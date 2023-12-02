from fractions import Fraction

import pytest

import gmpy2
from gmpy2 import mpz, mpq, mpfr, mpc, get_context, square, add, digits, xmpz
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


def test_digits():
    z2 = mpz(5)

    pytest.raises(TypeError, lambda: digits())
    pytest.raises(TypeError, lambda: digits(5, 5, 4, 5))

    assert digits(z2) == '5'
    assert digits(z2, 2) == '101'

    pytest.raises(TypeError, lambda: digits(z2, 2, 5))

    assert digits(mpq(3,5)) == '3/5'
    assert digits(mpq(3,5), 4) == '3/11'
    assert digits(mpfr(3,5), 4) == ('300', 1, 5)
    assert digits(mpfr(3,5), 4, 5) == ('30000', 1, 5)
    assert digits(complex(5,5), 4, 5) == (('11000', 2, 53), ('11000', 2, 53))

    pytest.raises(TypeError, lambda: digits('string', 4, 5))


def test_abs():
    a = mpz(123)
    b = abs(a)

    assert a is b

    a = xmpz(123)
    b = abs(a)

    assert a is not b

    a = mpz(-123)
    b = abs(a)

    assert b == mpz(123)
    assert a is not b
    assert a == mpz(-123)

    a = mpq(12,7)
    b = abs(a)

    assert a is b

    a = mpq(-12,7)
    b = abs(a)

    assert b == mpq(12,7)
    assert a == mpq(-12,7)

    a = xmpz(-123)
    b = abs(a)

    assert a == xmpz(123)
    assert b is None

    b = abs(a)

    assert b is None

    a = mpfr(1.0)
    b = abs(a)

    assert a is not b
    assert abs(mpfr(1, precision=100)) == mpfr('1.0')

    ctx = gmpy2.get_context()
    ctx.clear_flags()

    assert gmpy2.is_nan(abs(mpfr('nan')))
    assert ctx.invalid

    ctx.clear_flags()

    assert abs(mpfr('inf')) == mpfr('inf')
    assert ctx.overflow is False

    ctx.clear_flags()

    assert abs(mpfr('-inf')) == mpfr('inf')
    assert ctx.overflow is False

    assert abs(mpc(-1,0)) == mpfr('1.0')
    assert abs(-1+0j) == 1.0
    assert abs(mpc(1,1)) == mpfr('1.4142135623730951')

    ctx = gmpy2.get_context()
    ctx.clear_flags()

    c = mpc('nan+0j')

    assert gmpy2.is_nan(c.real) and c.imag == 0.0
    assert ctx.invalid
    assert gmpy2.is_nan(abs(c))
    assert ctx.invalid

    ctx.clear_flags()

    assert gmpy2.is_nan(abs(mpc('nanj')))
    assert ctx.invalid

    ctx.clear_flags()

    assert abs(mpc('inf+10j')) == mpfr('inf')
    assert ctx.overflow is False

    ctx.clear_flags()

    assert abs(mpc('-infj')) == mpfr('inf')
    assert ctx.overflow is False

    a = mpc('nan+infj')
    ctx.clear_flags()

    assert abs(a) == mpfr('inf')
    assert ctx.overflow is False
    assert ctx.invalid is False

    a=mpc('-inf+nanj')
    ctx.clear_flags()

    assert abs(a) == mpfr('inf')
    assert ctx.overflow is False
    assert ctx.invalid is False

    ctx = gmpy2.context()

    pytest.raises(TypeError, lambda: ctx.abs('a'))

    assert ctx.abs(-1) == mpz(1)
    assert ctx.abs(0) == mpz(0)
    assert ctx.abs(1) == mpz(1)
    assert ctx.abs(mpz(8)) == mpz(8)
    assert ctx.abs(mpz(-8)) == mpz(8)
    assert ctx.abs(-1.0) == mpfr('1.0')
    assert ctx.abs(mpfr(-2)) == mpfr('2.0')
    assert ctx.abs(2+3j) == mpfr('3.6055512754639891')
    assert ctx.abs(mpc(2+3j)) == mpfr('3.6055512754639891')
    assert ctx.abs(Fraction(1,2)) == mpq(1,2)
    assert ctx.abs(Fraction(-1,2)) == mpq(1,2)


def test_muldiv_2exp():
    ctx = gmpy2.get_context()
    r = mpfr(7.6)
    z = mpz(3)
    c = mpc(4,4)

    assert gmpy2.mul_2exp(r, z) == mpfr('60.799999999999997')
    assert gmpy2.mul_2exp(r, 3) == mpfr('60.799999999999997')

    assert gmpy2.mul_2exp(r, -5) == mpfr('0.23749999999999999')
    assert gmpy2.mul_2exp(1, -5) == mpfr('0.03125')  # issue 328

    pytest.raises(OverflowError, lambda: gmpy2.mul_2exp(z, 10**100))
    pytest.raises(OverflowError, lambda: gmpy2.mul_2exp(z, -10**100))
    pytest.raises(TypeError, lambda: gmpy2.mul_2exp(z, r))
    pytest.raises(TypeError, lambda: gmpy2.mul_2exp('not', 5))
    pytest.raises(TypeError, lambda: ctx.mul_2exp(r, z, 45))

    assert ctx.mul_2exp(c, z) == mpc('32.0+32.0j')
    assert ctx.mul_2exp(c, -5) == mpc('0.125+0.125j')

    pytest.raises(OverflowError, lambda: ctx.mul_2exp(c, 10**100))
    pytest.raises(OverflowError, lambda: ctx.mul_2exp(c, -10**100))

    assert ctx.mul_2exp(r, 0) == mpfr('7.5999999999999996')

    assert gmpy2.div_2exp(r, z) == mpfr('0.94999999999999996')
    assert gmpy2.div_2exp(r, 3) == mpfr('0.94999999999999996')
    assert gmpy2.div_2exp(r, -5) == mpfr('243.19999999999999')

    pytest.raises(OverflowError, lambda: gmpy2.div_2exp(r, 10**100))
    pytest.raises(OverflowError, lambda: gmpy2.div_2exp(r, -10**100))
    pytest.raises(TypeError, lambda: gmpy2.div_2exp(z, r))
    pytest.raises(TypeError, lambda: gmpy2.div_2exp('not', 5))
    pytest.raises(TypeError, lambda: ctx.div_2exp(r, z, 45))

    assert ctx.div_2exp(c, z) == mpc('0.5+0.5j')
    assert ctx.div_2exp(c, -5) == mpc('128.0+128.0j')

    pytest.raises(OverflowError, lambda: ctx.div_2exp(c, 10**100))
    pytest.raises(OverflowError, lambda: ctx.div_2exp(c, -10**100))
