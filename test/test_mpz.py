from hypothesis import given, example, settings
from hypothesis.strategies import booleans, integers, sampled_from
from pytest import raises

from gmpy2 import mpz


def test_mpz_to_bytes_interface():
    x = mpz(1)
    with raises(TypeError):
        x.to_bytes(1, 2, 3)
    with raises(TypeError):
        x.to_bytes(1, 2)
    with raises(TypeError):
        x.to_bytes('spam')
    with raises(TypeError):
        x.to_bytes(a=1, b=2, c=3, d=4)
    with raises(TypeError):
        x.to_bytes(2, length=2)
    with raises(TypeError):
        x.to_bytes(2, 'big', byteorder='big')
    with raises(TypeError):
        x.to_bytes(spam=1)

    with raises(ValueError):
        x.to_bytes(2, 'spam')
    with raises(ValueError):
        x.to_bytes(-1)

    assert x.to_bytes(2) == x.to_bytes(length=2)
    assert x.to_bytes(2, byteorder='little') == x.to_bytes(2, 'little')

    assert x.to_bytes() == x.to_bytes(1)
    assert x.to_bytes() == x.to_bytes(1, 'big')
    assert x.to_bytes() == x.to_bytes(signed=False)


@settings(max_examples=1000)
@given(integers(), integers(min_value=0, max_value=10000),
       sampled_from(['big', 'little']), booleans())
@example(0, 0, 'big', False)
@example(0, 0, 'little', False)
@example(0, 1, 'big', False)
@example(128, 1, 'big', True)
@example(128, 1, 'little', True)
@example(-129, 1, 'big', True)
@example(-129, 1, 'little', True)
@example(-1, 0, 'big', True)
@example(-2, 0, 'big', True)
@example(-2, 0, 'little', True)
@example(42, 1, 'big', False)
@example(42, 1, 'little', False)
@example(42, 3, 'big', False)
@example(42, 3, 'little', False)
@example(1000, 2, 'big', False)
@example(1000, 4, 'big', False)
@example(-2049, 1, 'big', True)
def test_to_bytes_mpz_vs_int(x, length, byteorder, signed):
    try:
        rx = x.to_bytes(length, byteorder, signed=signed)
    except OverflowError:
        with raises(OverflowError):
            mpz(x).to_bytes(length, byteorder, signed=signed)
    else:
        assert rx == mpz(x).to_bytes(length, byteorder, signed=signed)
