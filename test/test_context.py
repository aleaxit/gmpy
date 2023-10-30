import gmpy2
from gmpy2 import mpc, mpfr, mpz


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
