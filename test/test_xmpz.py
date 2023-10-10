from ctypes import memmove

import pytest

import gmpy2
from gmpy2 import mpz, xmpz, mpfr


def test_xmpz_digits():
    x = xmpz(6)
    assert x.digits() == '6'
    assert x.digits(2) == '0b110'
    assert x.digits(8) == '0o6'
    assert x.digits(16) == '0x6'
    assert xmpz(30).digits(16) == '0x1e'

    pytest.raises(ValueError, lambda: x.digits(0))
    pytest.raises(ValueError, lambda: x.digits(63))


def test_xmpz_limbs():
    x, y = xmpz(123456789), xmpz(0);
    x_limbs, num_limbs = x.limbs_read(), x.num_limbs()

    assert x_limbs != 0
    assert num_limbs > 0

    y_limbs = y.limbs_write(num_limbs)

    assert y_limbs != 0
    assert memmove(y_limbs, x_limbs, num_limbs * xmpz.limb_size) > 0

    y.limbs_finish(num_limbs)

    assert int(y) == 123456789

    x, y = xmpz(987654321), xmpz(0)
    x_limbs, num_limbs = x.limbs_read(), x.num_limbs()

    assert x_limbs != 0
    assert num_limbs > 0

    y_limbs = y.limbs_modify(num_limbs);

    assert y_limbs != 0
    assert memmove(y_limbs, x_limbs, num_limbs * xmpz.limb_size) > 0

    y.limbs_finish(num_limbs)

    assert int(y) == 987654321


def test_xmpz_attributes():
    x = xmpz(10)

    assert x.numerator == xmpz(10)
    assert x.denominator == xmpz(1)
    assert x.real == xmpz(10)


def test_xmpz_misc():
    z = xmpz(5)

    assert gmpy2.xbit_mask(7) == xmpz(127)
    assert gmpy2.xbit_mask(z) == xmpz(31)

    pytest.raises(TypeError, lambda: gmpy2.xbit_mask(4.5))

    assert not xmpz(0)
    assert z
    -z
    assert z == xmpz(-5)
    -z
    assert z == xmpz(5)
    +z
    assert z == xmpz(5)

    z2 = z.copy()

    assert z2 == z
    assert z2 is not z
    assert z.copy().make_mpz() == mpz(5)
    assert xmpz(-100).make_mpz() == mpz(-100)
    assert len(z) == 3
    assert len(xmpz(3000)) == 12


def test_xmpz_subscripts():
    x = xmpz(10)

    assert x[0] == x[2] == 0
    assert x[1] == x[3] == 1

    x[0] = 1

    assert x == xmpz(11)
    assert x[0:2] == mpz(3)
    assert x[0:] == mpz(11)
    assert x[:4] == mpz(11)

    x[:4] = 14

    assert x == xmpz(14)
    assert x[0] == 0

    pytest.raises(TypeError, lambda: x[mpfr('inf')])

    assert x[3:0] == mpz(0)
    assert x[-1] == 1

    x[-1] = 0

    assert x == xmpz(6)

    x[0:4] = 15

    assert x == xmpz(15)

    with pytest.raises(ValueError):
        x[0] = 16

    x[0:4] = 16

    assert x == xmpz(0)

    x[0:4] = 15

    assert x == xmpz(15)

    x[0:5] = 16

    assert x == xmpz(16)


def test_xmpz_iterators():
    x = xmpz(16)

    assert [b for b in x.iter_bits()] == [False, False, False, False, True]

    x = xmpz(30)

    assert [b for b in x.iter_bits()] == [False, True, True, True, True]
    assert [b for b in x.iter_bits(-1, 2)] == []
    assert [b for b in x.iter_bits(1, 3)] == [True, True]
    assert [b for b in x.iter_set()] == [1, 2, 3, 4]
    assert [b for b in x.iter_clear()] == [0]

    x = xmpz(10)

    assert [b for b in x.iter_clear()] == [0, 2]
