import pytest

import gmpy2
from gmpy2 import (context, get_context, ieee, local_context, mpc, mpfr, mpz,
                   set_context)


def test_context_abs():
    ctx = gmpy2.context()

    assert ctx.abs(-1) == mpz(1)
    assert ctx.abs(0) == mpz(0)
    assert ctx.abs(1) == mpz(1)
    assert ctx.abs(mpz(8)) == mpz(8)
    assert ctx.abs(mpz(-8)) == mpz(8)
    assert ctx.abs(-1.0) == mpfr('1.0')
    assert ctx.abs(mpfr(-2)) == mpfr('2.0')
    assert ctx.abs(2+3j) == mpfr('3.6055512754639891')
    assert ctx.abs(mpc(2+3j)) == mpfr('3.6055512754639891')


def test_context_ieee():
    ctx = ieee(32)

    assert (ctx.precision == 24 and ctx.emax == 128 and
            ctx.emin == -148 and ctx.subnormalize)

    ctx = ieee(64)

    assert (ctx.precision == 53 and ctx.emax == 1024 and
            ctx.emin == -1073 and ctx.subnormalize)

    ctx = ieee(128)

    assert (ctx.precision == 113 and ctx.emax == 16384 and
            ctx.emin == -16493 and ctx.subnormalize)

    ctx = ieee(256)

    assert (ctx.precision == 237 and ctx.emax == 262144 and
            ctx.emin == -262377 and ctx.subnormalize)

    pytest.raises(ValueError, lambda: ieee(-1))
    pytest.raises(TypeError, lambda: ieee("a"))

    set_context(ieee(32))

    assert gmpy2.const_pi().digits(2) == ('110010010000111111011011', 2, 24)

    set_context(ieee(64))

    assert gmpy2.const_pi().digits(2) == ('11001001000011111101101010100010001000010110100011000', 2, 53)

    set_context(ieee(128))

    assert gmpy2.const_pi().digits(2) == ('11001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111000', 2, 113)


def test_context():
    ctx = context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = context(precision=100)

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = context(real_prec=100)

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize and
            ctx.real_prec == 100)

    ctx = context(real_prec=100,imag_prec=200)

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize and
            ctx.real_prec == 100 and ctx.imag_prec == 200)


def test_get_context():
    set_context(context())

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()
    ctx.precision = 100

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    b = ctx.copy()
    b.precision = 200

    assert (b.precision == 200 and b.emax == 1073741823 and
            b.emin == -1073741823 and not b.subnormalize)
    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)


def test_local_context():
    set_context(context())

    with local_context() as ctx:
        assert ctx.precision == 53
        ctx.precision += 20
        assert ctx.precision == 73

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with local_context(ieee(64)) as ctx:
        assert (ctx.precision == 53 and ctx.emax == 1024 and
                ctx.emin == -1073 and ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with get_context() as ctx:
        assert ctx.precision == 53
        ctx.precision += 100
        assert ctx.precision == 153

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with local_context(precision=200) as ctx:
        assert ctx.precision == 200
        ctx.precision += 100
        assert ctx.precision == 300

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)
