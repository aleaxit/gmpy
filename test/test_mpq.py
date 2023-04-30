import decimal
import numbers
import pickle

from gmpy2 import mpq, mpz, cmp, cmp_abs
from supportclasses import q


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
    assert mpq.from_decimal(decimal.Decimal("5e-3")) == mpq(5, 1000)


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
