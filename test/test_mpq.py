import numbers

from gmpy2 import mpq, mpz


def test_mpq_as_integer_ratio():
    assert mpq(2, 3).as_integer_ratio() == (mpz(2), mpz(3))


def test_mpq_numbers_abc():
    assert isinstance(mpq(1, 3), numbers.Rational)
