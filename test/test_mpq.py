import numbers
import pickle
from decimal import Decimal

import pytest
from hypothesis import given, example, settings
from hypothesis.strategies import integers

from gmpy2 import mpq, mpz, cmp, cmp_abs, from_binary, to_binary
from supportclasses import a, b, c, d, q, z


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


def test_mpq_hash():
    hash(mpq(123456,1000)) == hash(Decimal('123.456'))
