from ctypes import memmove

import pytest

from gmpy2 import xmpz


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
