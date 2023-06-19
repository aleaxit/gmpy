import pytest

import gmpy2
from gmpy2 import (mpc, cmp, cmp_abs, nan, random_state, mpc_random,
                   to_binary, from_binary, get_context, is_nan)
from supportclasses import a, b, c, d


def test_mpc_cmp():
    pytest.raises(TypeError, lambda: cmp(mpc(1,2), mpc(3,4)))
    assert cmp_abs(mpc(1,2), mpc(3,4)) == -1
    assert cmp_abs(mpc(1,2), mpc(1,2)) == 0
    assert cmp_abs(mpc(3,4), mpc(1,2)) == 1
    gmpy2.get_context().clear_flags()
    assert gmpy2.get_context().erange is False
    assert cmp_abs(mpc(nan(),1), mpc(4.5)) == 0
    assert gmpy2.get_context().erange is True


def test_mpc_conversion():
    x = mpc(a)
    assert isinstance(x, mpc)
    assert x == 42+67j

    pytest.raises(TypeError, lambda: mpc(b))
    pytest.raises(TypeError, lambda: mpc(c))
    pytest.raises(TypeError, lambda: mpc(d))


def test_mpc_creation():
    ctx = gmpy2.get_context()
    ctx.clear_flags()
    a = mpc("1.2")
    assert a.rc == (-1, 0)
    assert ctx.inexact
    ctx.clear_flags()
    a = mpc("(1 2)")
    assert a == 1 + 2j
    assert a.rc == (0, 0)
    assert ctx.inexact is False
    ctx.clear_flags()
    a = mpc("1   + 2.1  j")
    assert a == 1 + 2.1j
    assert a.rc == (0, 1)
    assert ctx.inexact


def test_mpc_random():
    assert (mpc_random(random_state(42))
            == mpc('0.86555158787663011+0.4422082613292212j'))


def test_mpc_to_from_binary():
    x = mpc("0+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("1+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("-1+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0+1j")
    assert x == from_binary(to_binary(x))
    x = mpc("0-1j")
    assert x == from_binary(to_binary(x))
    x = mpc("inf+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0+infj")
    assert x == from_binary(to_binary(x))
    x = mpc("inf-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("inf+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("-inf+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("-inf-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("-inf+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("nan+0j")
    y = from_binary(to_binary(x))
    assert x.imag == y.imag
    assert is_nan(y.real)
    x = mpc("0+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("nan-infj")
    y = from_binary(to_binary(x))
    assert x.imag == y.imag
    assert is_nan(y.real)
    x = mpc("nan+nanj")
    y = from_binary(to_binary(x))
    assert is_nan(y.real)
    assert is_nan(y.imag)
    get_context().real_prec=100
    get_context().imag_prec=110
    assert (from_binary(to_binary(mpc("1.3-4.7j"))) ==
            mpc('1.2999999999999999999999999999994-4.7000000000000000000000000000000025j',
                (100,110)))
