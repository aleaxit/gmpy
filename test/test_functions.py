import pytest

import gmpy2
from gmpy2 import (can_round, fac, fma, fmma, fmms, fms, get_exp, ieee, is_nan,
                   maxnum, minnum, mpc, mpfr, mpq, mpz, powmod, powmod_sec,
                   root, rootn, set_exp, zero)


def test_root():
    assert root(zero(1), 2) == mpfr('0.0')
    assert root(zero(-1), 2) == mpfr('-0.0')
    assert root(zero(-1), 3) == mpfr('-0.0')
    assert root(zero(-1), 4) == mpfr('-0.0')
    assert root(zero(-1), 5) == mpfr('-0.0')
    assert root(zero(-1), 6) == mpfr('-0.0')
    assert root(2, 2) == mpfr('1.4142135623730951')
    assert root(mpz(2), 2) == mpfr('1.4142135623730951')
    assert root(mpq(2), 2) == mpfr('1.4142135623730951')
    assert root(mpfr(2), 2) == mpfr('1.4142135623730951')
    pytest.raises(TypeError, lambda: root(mpc(2), 2))
    assert is_nan(root(-2, 2))
    pytest.raises(OverflowError, lambda: root(2, -2))
    pytest.raises(TypeError, lambda: root(2, 0.5))


def test_rootn():
    assert rootn(zero(1), 2) == mpfr('0.0')
    assert rootn(zero(-1), 2) == mpfr('0.0')
    assert rootn(zero(-1), 3) == mpfr('-0.0')
    assert rootn(zero(-1), 4) == mpfr('0.0')
    assert rootn(zero(-1), 5) == mpfr('-0.0')
    assert rootn(zero(-1), 6) == mpfr('0.0')


def test_maxnum():
    a = mpfr("12.34")
    b = mpfr("45.67")
    nan = mpfr("nan")
    inf = mpfr("inf")

    assert maxnum(a, b) == mpfr('45.670000000000002')
    assert maxnum(b, a) == mpfr('45.670000000000002')
    assert maxnum(a, -b) == mpfr('12.34')
    assert maxnum(a, 123456) == mpfr('123456.0')
    assert maxnum(12345678901234567890, a) == mpfr('1.2345678901234567e+19')
    assert maxnum(0, -1) == mpfr('0.0')
    assert maxnum(1, inf) == mpfr('inf')
    assert maxnum(1, -inf) == mpfr('1.0')
    assert maxnum(nan, a) == mpfr('12.34')
    assert maxnum(a, nan) == mpfr('12.34')
    assert maxnum(nan, inf) == mpfr('inf')
    assert maxnum(nan, -inf) == mpfr('-inf')
    assert is_nan(maxnum(nan, nan))


def test_minnum():
    a = mpfr("12.34")
    b = mpfr("45.67")
    nan = mpfr("nan")
    inf = mpfr("inf")
    minf = mpfr("-inf")

    assert minnum(a, b) == mpfr('12.34')
    assert minnum(b, a) == mpfr('12.34')
    assert minnum(1, inf) == mpfr('1.0')
    assert minnum(minf, a) == mpfr('-inf')
    assert minnum(nan, inf) == mpfr('inf')
    assert is_nan(minnum(nan, nan))


def test_fused():
    assert fma(2,3,4) == mpz(10)
    assert fma(2,3,-4) == mpz(2)
    assert fma(2.0,3,-4) == mpfr('2.0')
    assert fma(2,3.0,-4) == mpfr('2.0')
    assert fma(2,3,-4.0) == mpfr('2.0')
    assert fma(2,mpfr(3),-4.0) == mpfr('2.0')
    assert fma(mpc(2),mpfr(3),-4.0) == mpc('2.0+0.0j')
    assert fms(2,3,4) == mpz(2)
    assert fms(2,3,-4) == mpz(10)

    assert ieee(128).fma(7,1/7,-1) == mpfr('-5.55111512312578270211815834045410156e-17',113)
    assert ieee(128).fma(7,mpq(1,7),-1) == mpq(0,1)

    pytest.raises(TypeError, lambda: fma(1,2,"r"))

    assert fma(1,2,mpq(3,4)) == mpq(11,4)
    assert fms(1,2,mpq(3,4)) == mpq(5,4)
    assert fms(1,mpfr(2),3) == mpfr('-1.0')
    assert fms(1,mpc(2),3) == mpc('-1.0+0.0j')

    assert fmma(2,3,4,5) == mpz(26)
    assert fmma(2,3,-4,5) == mpz(-14)
    assert fmma(2.0,3,-4, mpq(5)) == mpfr('-14.0')
    assert fmma(2,3.0,-4,5) == mpfr('-14.0')
    assert fmma(2,3,-4.0,5) == mpfr('-14.0')
    assert fmma(2,mpfr(3),-4.0,5) == mpfr('-14.0')

    pytest.raises(TypeError, lambda: fmma(mpc(2),mpfr(3),-4.0,5))

    assert fmms(2,3,4,5) == mpz(-14)
    assert fmms(2,3,-4,5) == mpz(26)
    assert fmms(2, 3, mpq(1, 2), 5) == mpq(7,2)
    assert fmms(2, 3, mpfr(1.2), 1) == mpfr('4.7999999999999998')

    assert ieee(128).fmma(7,1/7,-1,3/11) == mpfr('0.727272727272727237401994671017746441',113)
    assert ieee(128).fmma(7,mpq(1,7),-1,mpq(3,11)) == mpq(8,11)


def test_trigonometric():
    assert gmpy2.acos(mpc(0.2, 0.2)) == mpc('1.3735541886535356-0.20256635782456389j')
    assert gmpy2.acos(mpc(0.2, 0.2)) == gmpy2.acos(complex(0.2, 0.2))

    assert gmpy2.asin(mpc(0.2,0.2)) == mpc('0.1972421381413611+0.20256635782456389j')
    assert gmpy2.asin(mpc(2.0,0.2)) == mpc('1.4560834209500821+1.3245636864399635j')
    assert gmpy2.asin(mpc(0.2,0.2)) == gmpy2.asin(complex(0.2,0.2))

    assert gmpy2.atan(mpc(2.0, 2.0)) == mpc('1.311223269671635+0.23887786125685909j')
    assert gmpy2.atan(mpc(2.0, 2.0)) == gmpy2.atan(complex(2.0, 2.0))

    c = mpc(2,3)

    assert gmpy2.cos(c) == mpc('-4.189625690968807-9.109227893755337j')

    assert gmpy2.sin(c) == mpc('9.1544991469114301-4.1689069599665647j')

    assert gmpy2.sin_cos(c) == (mpc('9.1544991469114301-4.1689069599665647j'), mpc('-4.189625690968807-9.109227893755337j'))
    assert gmpy2.sin_cos(c) == gmpy2.sin_cos(complex(2,3))
    assert gmpy2.sin_cos(c) == (gmpy2.sin(c), gmpy2.cos(c))

    assert gmpy2.tan(mpc(4,5)) == mpc('8.9834776469715613e-05+1.0000132074347847j')

    assert gmpy2.atanh(mpc(2.0, 3.0)) == mpc('0.14694666622552977+1.3389725222944935j')
    assert gmpy2.atanh(mpc(2.0, 3.0)) == gmpy2.atanh(complex(2, 3))

    assert gmpy2.tanh(mpc(4,5)) == mpc('1.0005630461157933-0.00036520305451130409j')


def test_get_exp():
    ctx = gmpy2.get_context()
    ctx.trap_erange = True

    pytest.raises(gmpy2.RangeError, lambda: get_exp(mpfr('inf')))


def test_set_exp():
    pytest.raises(ValueError, lambda: set_exp(mpfr('1.0'), int(fac(100))))

    gmpy2.set_context(gmpy2.ieee(32))
    ctx = gmpy2.get_context()
    ctx.trap_erange = True

    pytest.raises(gmpy2.RangeError, lambda: set_exp(mpfr('1.0'), 1000))

    ctx.trap_erange = False
    assert set_exp(mpfr('1.0'), 1000) == mpfr('1.0')


def test_can_round():
    pytest.raises(TypeError, lambda: can_round(mpfr('1.1'), 10, "spam"))
    pytest.raises(ValueError, lambda: can_round(mpfr('1.1'), 10, 111, 111, 111))
    pytest.raises(ValueError, lambda: can_round(mpfr('1.1'), 10, 1, 111, 111))
    pytest.raises(ValueError, lambda: can_round(mpfr('1.1'), 10, 1, 1, -111))

    x = mpfr('-1.112')

    assert can_round(x, 10, 1, 1, 1)
    assert not can_round(x, 10, 1, 1, 10)


def test_powmod():
    z1, z2 = mpz(5), mpz(2)
    q = mpq(2,3)

    assert powmod(z1, z2, 4) == pow(z1, z2, 4)

    pytest.raises(TypeError, lambda: powmod(z1))
    pytest.raises(TypeError, lambda: powmod(z1, q, 4))


def test_powmod_sec():
    assert powmod_sec(3,3,7) == mpz(6)
    assert powmod_sec(-3,3,7) == mpz(1)
    assert powmod(-3,3,7) == mpz(1)
    assert powmod(3,-3,7) == mpz(6)

    pytest.raises(ValueError, lambda: powmod_sec(3,-3,7))
    pytest.raises(ValueError, lambda: powmod_sec(3,4,8))
