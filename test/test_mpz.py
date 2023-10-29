import math
import numbers
import pickle

from hypothesis import assume, example, given, settings
from hypothesis.strategies import booleans, integers, sampled_from
from pytest import mark, raises
from supportclasses import a, b, c, d, q, z

import gmpy2
from gmpy2 import (cmp, cmp_abs, from_binary, mp_version, mpc, mpfr, mpq, mpz,
                   mpz_random, mpz_rrandomb, mpz_urandomb, pack, random_state,
                   to_binary, unpack, xmpz)


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
@example(321, -123)
def test_mpz_arithmetics(i, z):
    i, z = map(mpz, [i, z])

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
        assert i % int(z) == i % z
        assert int(i) % z == i % z
        assert divmod(int(i), int(z)) == divmod(i, z)
        assert divmod(i, int(z)) == divmod(i, z)
        assert divmod(int(i), z) == divmod(i, z)

    if i:
        assert int(z) // int(i) == z // i
        assert int(z) % int(i) == z % i
        assert z % int(i) == z % i
        assert int(z) % i == z % i
        assert divmod(int(z), int(i)) == divmod(z, i)
        assert divmod(z, int(i)) == divmod(z, i)
        assert divmod(int(z), i) == divmod(z, i)


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

    assert float(mpz(1)) == 1.0
    raises(OverflowError, lambda: float(mpz(99**199)))
    assert mpz(xmpz(1)) == mpz(1)


def test_mpz_create():
    assert mpz() == mpz(0)
    assert mpz(0) == mpz(0)
    assert mpz(1) == mpz(1)
    assert mpz(-1) == mpz(-1)
    assert mpz(2**15-2) == mpz(32766)
    assert mpz(2**15-1) == mpz(32767)
    assert mpz(2**15) == mpz(32768)
    assert mpz(2**15+1) == mpz(32769)
    assert mpz(2**15+2) == mpz(32770)
    assert mpz(2**30-2) == mpz(1073741822)
    assert mpz(2**30-1) == mpz(1073741823)
    assert mpz(2**30) == mpz(1073741824)
    assert mpz(2**30+1) == mpz(1073741825)
    assert mpz(2**30+2) == mpz(1073741826)
    assert mpz(2**16-2) == mpz(65534)
    assert mpz(2**16-1) == mpz(65535)
    assert mpz(2**16) == mpz(65536)
    assert mpz(2**16+1) == mpz(65537)
    assert mpz(2**16+2) == mpz(65538)
    assert mpz(1000000000000) == mpz(1000000000000)
    assert mpz(-1000000000000) == mpz(-1000000000000)

    raises(ValueError, lambda: mpz(''))
    raises(ValueError, lambda: mpz('a'))

    assert mpz('a',16) == mpz(10)

    raises(ValueError, lambda: mpz('z',16))

    assert mpz('0b1101') == mpz(13)
    assert mpz('0b1101',2) == mpz(13)
    assert mpz('1101',2) == mpz(13)
    assert mpz('0b0010') == mpz(2)
    assert mpz('0b0010',2) == mpz(2)

    raises(ValueError, lambda: mpz('0b0b10',2))
    raises(ValueError, lambda: mpz('0b0b10'))
    raises(ValueError, lambda: mpz('0b0012'))

    assert mpz('0o0012') == mpz(10)
    assert mpz('0o0012',8) == mpz(10)
    assert mpz('12',8) == mpz(10)
    assert mpz('0x12') == mpz(18)
    assert mpz('0x12',16) == mpz(18)
    assert mpz('12',16) == mpz(18)
    assert mpz('-1') == mpz(-1)
    assert mpz('+1') == mpz(1)
    assert mpz('  0xA', base=0) == mpz(10)

    raises(ValueError, lambda: mpz(float('nan')))
    raises(OverflowError, lambda: mpz(float('inf')))
    raises(OverflowError, lambda: mpz(float('-inf')))
    raises(TypeError, lambda: mpz(12, base=16))

    assert mpz('12', base=16) == mpz(18)

    raises(ValueError, lambda: mpz('\xff'))
    raises(ValueError, lambda: mpz('\x0cf'))
    raises(ValueError, lambda: mpz('\0xff'))

    assert mpz(b'12') == mpz(12)

    raises(TypeError, lambda: mpz(None))
    raises(TypeError, lambda: mpz(None,base=10))
    raises(ValueError, lambda: mpz('99',base=100))
    raises(TypeError, lambda: mpz('99',base='a'))

    assert mpz('99',base=10) == mpz(99)
    assert mpz(xmpz(5)) == mpz(5)

    raises(ValueError, lambda: mpz('ы'))
    raises(ValueError, lambda: mpz(bytes('ы', encoding='utf-8')))


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


@settings(max_examples=1000)
@given(integers())
@example(0)
@example(1)
@example(-1)
@example(-2)
@example(123)
def test_mpz_hash(n):
    assert hash(mpz(n)) == hash(n)


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


def test_mpz_format():
    z1, z2 = mpz(-3), mpz(5)

    assert '{:<5}'.format(z1) == '-3   '
    assert '{:>+5}'.format(z2) == '   +5'

    raises(ValueError, lambda: '{:5+}'.format(z1))

    assert '{:>-4}'.format(z2) == '   5'
    assert '{:<-4}'.format(z1) == '-3  '

    raises(ValueError, lambda: '{:>4-}'.format(z1))
    raises(ValueError, lambda: '{:<4 }'.format(z1))

    assert '{:#x}'.format(z1) == '-0x3'
    assert '{:#o}'.format(z1) == '-0o3'

    raises(ValueError, lambda: '{:>5#}'.format(z1))
    raises(ValueError, lambda: '{:~}'.format(z1))


def test_mpz_digits():
    z1, z2 = mpz(-3), mpz(15)

    assert z1.digits() == '-3'
    assert z1.digits(2) == '-11'
    assert z1.digits(8) == '-3'
    assert z2.digits(16) == 'f'

    raises(ValueError, lambda: z1.digits(0))
    raises(ValueError, lambda: z1.digits(1))


def test_mpz_sub():
    a, b = mpz(123), mpz(456)
    c = 12345678901234567890

    assert a-1 == mpz(122)
    assert a-(-1) == mpz(124)
    assert 1-a == mpz(-122)
    assert (-1)-a == mpz(-124)
    assert a-c == -12345678901234567767
    assert c-a == 12345678901234567767
    assert a-(-c) == 12345678901234568013
    assert (-c)-a == -12345678901234568013
    assert a-b == mpz(-333)
    assert b-a == mpz(333)
    assert a-(-b) == mpz(579)
    assert (-b)-a == mpz(-579)
    assert a-z == mpz(121)
    assert gmpy2.sub(2,1) == mpz(1)

    ctx = gmpy2.context()

    assert ctx.sub(a,b) == a-b
    assert ctx.sub(c,c) == c-c
    assert ctx.sub(1, 1) == mpz(0)
    assert ctx.sub(a, 1) == mpz(122)
    assert ctx.sub(1, a) == mpz(-122)
    assert ctx.sub(a, mpq(0)) == mpq(123,1)
    assert ctx.sub(a, mpfr(0)) == mpfr('123.0')
    assert ctx.sub(a, mpc(0)) == mpc('123.0+0.0j')

    raises(TypeError, lambda: ctx.sub(1))
    raises(TypeError, lambda: ctx.sub(1,2,3))
    raises(TypeError, lambda: a-'b')
    raises(TypeError, lambda: 'b'-a)


def test_mpz_mul():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890

    assert mpz(2)*z == mpz(4)
    assert gmpy2.mul(2,1) == mpz(2)

    ctx = gmpy2.context()

    assert ctx.mul(a,b) == a*b
    assert ctx.mul(c,c) == c*c
    assert ctx.mul(a, mpq(1)) == mpq(123,1)
    assert ctx.mul(a, mpfr(1)) == mpfr('123.0')
    assert ctx.mul(a, mpc(1)) == mpc('123.0+0.0j')

    raises(TypeError, lambda: ctx.mul(1))
    raises(TypeError, lambda: ctx.mul(1,2,3))


def test_mpz_divmod():
    a = mpz(123)
    b = mpz(456)

    raises(TypeError, lambda: divmod(mpz(123),'a'))

    ctx = gmpy2.context()

    raises(TypeError, lambda: ctx.divmod('a',456))
    raises(TypeError, lambda: ctx.divmod(1,2,3))
    raises(ZeroDivisionError, lambda: ctx.divmod(456, 0))
    raises(TypeError, lambda: ctx.divmod(a,mpc(456)))
    raises(TypeError, lambda: divmod(mpz(1), mpc(1,2)))

    assert ctx.divmod(a,b) == (mpz(0), mpz(123))
    assert ctx.divmod(123,456) == (mpz(0), mpz(123))
    assert divmod(mpz(3), z) == (mpz(1), mpz(1))
    assert divmod(z, mpz(3)) == (mpz(0), mpz(2))


def test_mpz_pow():
    z1, z2 = mpz(5), mpz(2)
    ctx = gmpy2.get_context()

    assert z1 ** z2 == mpz(25)
    assert ctx.pow(z1, z2) == mpz(25)
    assert z1 ** -z2 == mpfr('0.040000000000000001')
    assert z1 ** 0 == mpz(1)
    assert mpz(0) ** 32 == mpz(0)
    assert mpz(-1) ** 32 == mpz(1)
    assert mpz(1) ** mpz(324) == mpz(1)
    assert mpz(0) ** 0 == mpz(1)
    assert mpz(-1) ** 3 == mpz(-1)
    assert z1 ** 2 == pow(z1, 2)
    assert pow(z1, 2, 19) == mpz(6)
    assert pow(z1, -2, 19) == mpz(16)

    raises(ValueError, lambda: pow(mpz(0), -2, 19))

    assert pow(z1, 2, -19) == mpz(-13)

    raises(ValueError, lambda: pow(5, 2, 0))
    raises(TypeError, lambda: ctx.pow(z1, 'invalid'))


def test_lucasu():
    assert gmpy2.lucasu(2,4,1) == mpz(1)

    raises(ValueError, lambda: gmpy2.lucasu(2,1,1))

    assert gmpy2.lucasu(2,4,8) == mpz(128)

    raises(TypeError, lambda: gmpy2.lucasu('a',4,8))
    raises(TypeError, lambda: gmpy2.lucasu(2,'b',8))
    raises(TypeError, lambda: gmpy2.lucasu(2,4,None))
    raises(ValueError, lambda: gmpy2.lucasu(mpz(2), mpz(1), mpz(7)))

    assert gmpy2.lucasu_mod(3,2,5,7) == mpz(3)
    assert gmpy2.lucasu_mod(3,2,555,777777777) == mpz(387104641)

    raises(ValueError, lambda: gmpy2.lucasu_mod(2,1,555,777777777))

    raises(ValueError, lambda: gmpy2.lucasv(2,1,4))
    raises(TypeError, lambda: gmpy2.lucasv('a',1,2))
    raises(TypeError, lambda: gmpy2.lucasv(4,'b',2))
    raises(TypeError, lambda: gmpy2.lucasv(4,3,'c'))

    assert gmpy2.lucasv(4,3,7) == mpz(2188)
    assert gmpy2.lucasv(4,3,8) == mpz(6562)

    assert gmpy2.lucasv_mod(4,3,55,123456) == mpz(35788)
    assert gmpy2.lucasv_mod(4,3,56,123456) == mpz(107362)
    assert gmpy2.lucasv_mod(4,3,57,123456) == mpz(75172)
