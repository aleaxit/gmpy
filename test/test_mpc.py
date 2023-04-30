import pytest

import gmpy2
from gmpy2 import mpc, cmp, cmp_abs, nan
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
