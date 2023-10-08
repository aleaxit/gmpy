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
