import math
import numbers
import pickle
from decimal import Decimal

from hypothesis import assume, given, example, settings
from hypothesis.strategies import booleans, integers, sampled_from
from pytest import raises, mark

from gmpy2 import (mpz, pack, unpack, cmp, cmp_abs, to_binary, from_binary,
                   random_state, mpz_random, mpz_urandomb, mpz_rrandomb,
                   mp_version)
from supportclasses import a, b, c, d, z, q


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
@example(-65281, 3, 'big', True)
@example(-65281, 3, 'little', True)
def test_mpz_to_bytes(x, length, byteorder, signed):
    try:
        rx = x.to_bytes(length, byteorder, signed=signed)
    except OverflowError:
        with raises(OverflowError):
            mpz(x).to_bytes(length, byteorder, signed=signed)
    else:
        assert rx == mpz(x).to_bytes(length, byteorder, signed=signed)


def test_mpz_from_bytes_interface():
    with raises(TypeError):
        mpz.from_bytes()
    with raises(TypeError):
        mpz.from_bytes(1, 2, 3)
    with raises(TypeError):
        mpz.from_bytes(b'', 2)
    with raises(TypeError):
        mpz.from_bytes(1)
    with raises(TypeError):
        mpz.from_bytes(b'', bytes=b'')
    with raises(TypeError):
        mpz.from_bytes(b'', 'big', byteorder='big')
    with raises(TypeError):
        mpz.from_bytes(b'', spam=1)
    with raises(TypeError):
        mpz.from_bytes(a=1, b=2, c=3, d=4)

    with raises(ValueError):
        mpz.from_bytes(b'', 'spam')

    assert mpz.from_bytes(b'\x01', byteorder='little') == mpz.from_bytes(b'\x01', 'little')

    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(bytes=b'\x01')
    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(b'\x01', 'big')
    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(b'\x01', signed=False)


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
@example(-1, 1, 'big', True)
@example(-1, 1, 'little', True)
@example(-2, 0, 'big', True)
@example(-2, 0, 'little', True)
@example(-1, 3, 'big', True)
@example(-2, 3, 'big', True)
@example(-2, 5, 'little', True)
def test_mpz_from_bytes(x, length, byteorder, signed):
    try:
        bytes = x.to_bytes(length, byteorder, signed=signed)
    except OverflowError:
        assume(False)
    else:
        rx = int.from_bytes(bytes, byteorder, signed=signed)
        assert rx == mpz.from_bytes(bytes, byteorder, signed=signed)
        assert rx == mpz.from_bytes(bytearray(bytes), byteorder, signed=signed)
        assert rx == mpz.from_bytes(list(bytes), byteorder, signed=signed)


def test_mpz_as_integer_ratio():
    assert mpz(3).as_integer_ratio() == (mpz(3), mpz(1))


def test_mpz_numbers_abc():
    assert isinstance(mpz(2), numbers.Integral)


def test_mpz_pickling():
    for proto in range(pickle.HIGHEST_PROTOCOL + 1):
        for x in [mpz(12346789), mpz(-12346789), mpz(0)]:
            assert pickle.loads(pickle.dumps(x, protocol=proto)) == x


@settings(max_examples=1000)
@given(integers(), integers())
@example(0, 0)
@example(0, 1)
@example(-11, 75)
@example(14, 105)
@example(64, 123456789012345678901234567890)
def test_mpz_arithmetics(i, z):
    assert int(i) + int(z) == i + z
    assert int(z) + int(i) == z + i

    assert int(i) - int(z) == i - z
    assert int(z) - int(i) == z - i

    assert int(i) * int(z) == i * z
    assert int(z) * int(i) == z * i

    # Test all permutations of floor division
    if z:
        assert int(i) // int(z) == i // z
        assert int(i) % int(z) == i % z
        assert divmod(int(i), int(z)) == divmod(i, z)

    if i:
        assert int(z) // int(i) == z // i
        assert int(z) % int(i) == z % i
        assert divmod(int(z), int(i)) == divmod(z, i)

@settings(max_examples=1000)
@given(integers(min_value=0),
       integers(min_value=1, max_value=100000))
def test_mpz_pack_unpack_bulk(x, n):
    lst = unpack(x, n)
    assert pack(lst, n) == x


def test_mpz_pack_unpack():
    x = mpz(0)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(2)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(3141592635)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1234567891234567890000000000000000000000000000000000000123)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1) << 500
    assert all((x == pack(unpack(x,i),i) for i in range(1,200)))
    x -= 1
    assert all((x == pack(unpack(x,i),i) for i in range(1,200)))

    raises(TypeError, lambda: pack(x))
    raises(TypeError, lambda: pack(1, 1))
    raises(TypeError, lambda: pack([mpz(1), mpz(-1)], 2))

    raises(TypeError, lambda: unpack(x))
    raises(TypeError, lambda: unpack([], 1))
    raises(ValueError, lambda: unpack(-1, 1))


def test_mpz_cmp():
    assert cmp(0, mpz(0)) == 0
    assert cmp(1, mpz(0)) == 1
    assert cmp(0, mpz(1)) == -1
    assert cmp(-1, mpz(0)) == -1
    assert cmp(0, mpz(-1)) == 1

    assert cmp_abs(mpz(0), 0) == 0
    assert cmp_abs(mpz(1), 0) == 1
    assert cmp_abs(mpz(0), 1) == -1
    assert cmp_abs(mpz(-1), 0) == 1
    assert cmp_abs(mpz(0), -1) == -1

    a = mpz(-10)
    assert cmp_abs(a, 0) == 1
    assert a == mpz(-10)
    assert cmp_abs(100, a) == 1
    assert a == mpz(-10)

    assert cmp(mpz(2), z) == 0
    assert cmp(z, mpz(3)) == -1
    assert cmp(mpz(1), q) == -1
    assert cmp(mpz(1), mpz(q)) == 0


def test_mpz_conversion():
    x = mpz(a)
    assert isinstance(x, mpz)
    assert x == 1
    raises(TypeError, lambda: mpz(b))
    raises(TypeError, lambda: mpz(c))
    raises(TypeError, lambda: mpz(d))


@given(integers())
@example(0)
def test_mpz_conversion_bulk(n):
    assert int(mpz(n)) == n


@settings(max_examples=1000)
@given(integers())
@example(0)
@example(1)
@example(-1)
@example(123456789123456789)
def test_mpz_to_from_binary(n):
    x = mpz(n)
    assert x == from_binary(to_binary(x))


def test_mpz_hash():
    assert hash(mpz(123)) == hash(Decimal(123))


def test_mpz_ceil():
    a = mpz(123)
    assert math.ceil(a) == a
    assert math.ceil(a) is a


def test_mpz_floor():
    a = mpz(123)
    assert math.floor(a) == a
    assert math.floor(a) is a


def test_mpz_trunc():
    a = mpz(123)
    assert math.trunc(a) == a
    assert math.trunc(a) is a


def test_mpz_round():
    assert round(mpz(123456), 2) == mpz(123456)
    assert round(mpz(123456), -22) == mpz(0)
    assert round(mpz(123456), -2) == mpz(123500)
    assert round(mpz(123456), -1) == mpz(123460)
    assert round(mpz(123455), -1) == mpz(123460)
    assert round(mpz(123454), -1) == mpz(123450)
    assert round(mpz(123445), -1) == mpz(123440)
    assert round(mpz(123445)) == mpz(123445)


def test_mpz_random():
    r1 = random_state(42)
    r2 = random_state(42)

    assert mpz_random(r1, 2**88) == mpz(171378365038768291737094841)
    assert mpz_random(r2, 2**88) == mpz(171378365038768291737094841)
    assert mpz_random(r1, 2**88) == mpz(62749575961297098445301393)
    assert mpz_random(r2, 2**88) == mpz(62749575961297098445301393)


def test_mpz_urandomb():
    assert (mpz_urandomb(random_state(42), 64).digits(2) ==
            '1100100011011011101100101001100100111110010111011100101010111001')


def test_mpz_rrandomb():
    assert (mpz_rrandomb(random_state(42), 64).digits(2) ==
            '1111111111111111111111111100000000111111111111111111000000000000')

@mark.skipif(mp_version() < "GMP 6.3.0", reason="requires GMP 6.3.0 or higher")
def test_prev_prime():
    # Imported here as symbol won't exist if mp_version() < 6.3.0
    from gmpy2 import prev_prime
    assert prev_prime(3) == mpz(2)
    assert prev_prime(4) == mpz(3)
    assert prev_prime(5) == mpz(3)
    assert prev_prime(6) == mpz(5)
    assert prev_prime(1000004) == mpz(1000003)
    assert prev_prime(1000033) == mpz(1000003)

    with raises(ValueError):
        prev_prime(-100)

    with raises(ValueError):
        prev_prime(2)

    with raises(TypeError):
        prev_prime('a')

    with raises(TypeError):
        prev_prime(4.5)
