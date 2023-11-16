import pytest

import gmpy2
from gmpy2 import (can_round, check_range, copy_sign, f2q, fac, fma, fmma,
                   fmms, fms, from_binary, get_emax_max, get_emin_min, get_exp,
                   ieee, inf, is_bpsw_prp, is_euler_prp,
                   is_extra_strong_lucas_prp, is_fermat_prp, is_fibonacci_prp,
                   is_finite, is_infinite, is_lucas_prp, is_nan,
                   is_selfridge_prp, is_strong_bpsw_prp, is_strong_lucas_prp,
                   is_strong_prp, is_strong_selfridge_prp, is_zero, maxnum,
                   minnum, mpc, mpfr, mpfr_from_old_binary, mpq,
                   mpq_from_old_binary, mpz, mpz_from_old_binary, nan, norm,
                   phase, polar, powmod, powmod_sec, proj, rect, root,
                   root_of_unity, rootn, set_exp, set_sign, sign, zero)


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


def test_is_fermat_prp():
    assert is_fermat_prp(12345,2) is False
    assert is_fermat_prp(113,2)
    assert is_fermat_prp(1234,2) is False

    pytest.raises(TypeError, lambda: is_fermat_prp(1234,'a'))
    pytest.raises(TypeError, lambda: is_fermat_prp(1234, 2, 3))
    pytest.raises(ValueError, lambda: is_fermat_prp(113, 1))
    pytest.raises(ValueError, lambda: is_fermat_prp(-113, 3))
    pytest.raises(ValueError, lambda: is_fermat_prp(339, 3))

    assert is_fermat_prp(mpz(12345),2) is False
    assert is_fermat_prp(113,mpz(2))


def test_is_euler_prp():
    assert is_euler_prp(12345,2) is False
    assert is_euler_prp(113,2)
    assert is_euler_prp(1234,2) is False

    pytest.raises(TypeError, lambda: is_euler_prp(1234,'a'))
    pytest.raises(TypeError, lambda: is_euler_prp(1234, 2, 3))
    pytest.raises(ValueError, lambda: is_euler_prp(113, 1))
    pytest.raises(ValueError, lambda: is_euler_prp(-113, 3))
    pytest.raises(ValueError, lambda: is_euler_prp(339, 3))

    assert is_euler_prp(mpz(12345),2) is False
    assert is_euler_prp(113,mpz(2))


def test_is_strong_prp():
    assert is_strong_prp(12345,2) is False
    assert is_strong_prp(113,2)
    assert is_strong_prp(1234,2) is False

    pytest.raises(TypeError, lambda: is_strong_prp(1234,'a'))
    pytest.raises(TypeError, lambda: is_strong_prp(1234, 2, 3))
    pytest.raises(ValueError, lambda: is_strong_prp(113, 1))
    pytest.raises(ValueError, lambda: is_strong_prp(-113, 3))
    pytest.raises(ValueError, lambda: is_strong_prp(339, 3))

    assert is_strong_prp(mpz(12345),2) is False
    assert is_strong_prp(113,mpz(2))


def test_is_fibonacci_prp():
    assert is_fibonacci_prp(12345, 3, 1) is False
    assert is_fibonacci_prp(113, 3, 1)
    assert is_fibonacci_prp(12345, 3, -1) is False
    assert is_fibonacci_prp(113, 3, -1)

    pytest.raises(ValueError, lambda: is_fibonacci_prp(113, 3, 2))
    pytest.raises(TypeError, lambda: is_fibonacci_prp('a', 3, 2))
    pytest.raises(ValueError, lambda: is_fibonacci_prp(113, 2, 1))

    assert is_fibonacci_prp(113, 2, -1)


def test_is_lucas_prp():
    assert is_lucas_prp(12345, 5, 2) is False
    assert is_lucas_prp(113, 5, 2)

    pytest.raises(ValueError, lambda: is_lucas_prp(12345, 3, 5))


def test_is_is_stronglucas_prp():
    assert is_strong_lucas_prp(12345, 5, 2) is False
    assert is_strong_lucas_prp(113, 5, 2)

    pytest.raises(ValueError, lambda: is_strong_lucas_prp(12345, 3, 5))


def test_is_extra_strong_lucas_prp():
    assert is_extra_strong_lucas_prp(12345, 9) is False
    assert is_extra_strong_lucas_prp(113, 5)

    pytest.raises(ValueError, lambda: is_extra_strong_lucas_prp(12345, 3))


def test_is_selfridge_prp():
    assert is_selfridge_prp(12345) is False
    assert is_selfridge_prp(113)


def test_is_strong_selfridge_prp():
    assert is_strong_selfridge_prp(12345) is False
    assert is_strong_selfridge_prp(113)


def test_is_bpsw_prp():
    assert is_bpsw_prp(12345) is False
    assert is_bpsw_prp(113)


def test_is_strong_bpsw_prp():
    assert is_strong_bpsw_prp(12345) is False
    assert is_strong_bpsw_prp(113)


def test_mpz_from_old_binary():
    assert gmpy2.mpz_from_old_binary(b'\x15\xcd[\x07') == mpz(123456789)
    assert gmpy2.mpz_from_old_binary(b'\x15\xcd[\x07\xff') == mpz(-123456789)

    pytest.raises(TypeError, lambda: mpz_from_old_binary(1))


def test_mpq_from_old_binary():
    assert mpq_from_old_binary(b'\x01\x00\x00\x00)\x98') == mpq(41,152)
    assert mpq_from_old_binary(b'\x01\x00\x00\x80)\x98') == mpq(-41,152)

    pytest.raises(TypeError, lambda: mpq_from_old_binary(1))
    pytest.raises(ValueError, lambda: mpq_from_old_binary(b'aa'))
    pytest.raises(ValueError, lambda: mpq_from_old_binary(b'aaaaaaaaa'))


def test_mpfr_from_old_binary():
    assert mpfr_from_old_binary(b'\x085\x00\x00\x00\x02\x00\x00\x0009\xac\xcc\xcc\xcc\xcc\xcc\xd0') == mpfr('12345.674999999999')
    assert mpfr_from_old_binary(b'\t5\x00\x00\x00\x02\x00\x00\x0009\xac\xcc\xcc\xcc\xcc\xcc\xd0') == mpfr('-12345.674999999999')
    assert mpfr_from_old_binary(b'\n5\x00\x00\x00\x06\x00\x00\x00\x01\x14\xb3\x7fKQ\xf7\x0en') == mpfr('1.5e-17')

    pytest.raises(TypeError, lambda: mpfr_from_old_binary(1))
    pytest.raises(ValueError, lambda: mpfr_from_old_binary(b'aaaaa'))

    assert mpfr_from_old_binary(b'\x04') == mpfr('0.0')


def test_from_binary():
    pytest.raises(TypeError, lambda: from_binary(1))
    pytest.raises(ValueError, lambda: from_binary(b'a'))


def test_phase():
    pytest.raises(TypeError, lambda: phase())
    pytest.raises(TypeError, lambda: phase(3))

    assert phase(mpc(4,5)) == mpfr('0.89605538457134393')
    assert ieee(64).phase(mpc(4,5)) == mpfr('0.89605538457134393')


def test_root_of_unity():
    assert root_of_unity(1,1) == mpc('1.0+0.0j')
    assert root_of_unity(1,2) == mpc('1.0+0.0j')
    assert root_of_unity(2,1) == mpc('-1.0+0.0j')
    assert root_of_unity(3,1) == mpc('-0.5+0.8660254037844386j')
    assert root_of_unity(3,2) == mpc('-0.5-0.8660254037844386j')
    assert root_of_unity(3,3) == mpc('1.0+0.0j')
    assert ieee(128).root_of_unity(3,1) == mpc('-0.5+0.866025403784438646763723170752936161j',(113,113))

    pytest.raises(TypeError, lambda: ieee(128).root_of_unity())
    pytest.raises(TypeError, lambda: ieee(128).root_of_unity('a','b'))


def test_norm():
    pytest.raises(TypeError, lambda: norm())
    pytest.raises(TypeError, lambda: norm(2))

    assert norm(mpc(1,2)) == mpfr('5.0')
    assert ieee(32).norm(mpc(1,2)) == mpfr('5.0',24)


def test_polar():
    pytest.raises(TypeError, lambda: polar())
    pytest.raises(TypeError, lambda: polar(5))
    assert polar(mpc(1,1)) == (mpfr('1.4142135623730951'), mpfr('0.78539816339744828'))


def test_rect():
    pytest.raises(TypeError, lambda: rect())
    pytest.raises(TypeError, lambda: rect(1))

    assert rect(1,1) == mpc('0.54030230586813977+0.8414709848078965j')


def test_proj():
    pytest.raises(TypeError, lambda: proj())
    pytest.raises(TypeError, lambda: proj(1))

    assert proj(mpc(1,1)) == mpc('1.0+1.0j')
    assert proj(mpc(1,2)) == mpc('1.0+2.0j')


def test_is_zero():
    assert is_zero(mpc("0+0j"))
    assert is_zero(mpc("1+0j")) is False
    assert is_zero(mpc("1+1j")) is False
    assert is_zero(mpc("0+1j")) is False

def test_is_nan():
    assert is_nan(mpc("nan+1j"))
    assert is_nan(mpc("1+nanj"))
    assert is_nan(mpc("nan+nanj"))
    assert is_nan(mpc("1+1j")) is False

def test_is_infinite():
    assert is_infinite(mpc("inf+1j"))
    assert is_infinite(mpc("-inf+1j"))
    assert is_infinite(mpc("1+infj"))
    assert is_infinite(mpc("1-infj"))
    assert is_infinite(mpc("inf-infj"))
    assert is_infinite(mpc("1+1j")) is False

def test_is_finite():
    assert is_finite(mpc("0+0j"))
    assert is_finite(mpc("nan+0j")) is False
    assert is_finite(mpc("0+nanj")) is False
    assert is_finite(mpc("0+infj")) is False
    assert is_finite(mpc("inf+3j")) is False


def test_f2q():
    a = mpfr('123.456')

    pytest.raises(TypeError, lambda: f2q('a'))
    pytest.raises(TypeError, lambda: f2q(1,2,3,4))

    assert f2q(a,0.1) == mpz(123)
    assert f2q(a,0.01) == mpz(123)
    assert f2q(a,0.001) == mpq(247,2)
    assert f2q(a,0.0001) == mpq(1358,11)
    assert f2q(a,0.00001) == mpq(7037,57)
    assert f2q(a,0.000001) == mpq(15432,125)
    assert f2q(a,0.0000001) == mpq(15432,125)
    assert f2q(a) == mpq(15432,125)
    assert f2q(2.50000000000008) == mpq(15637498706148,6254999482459)
    assert f2q(2.5000000000000) == mpq(5,2)
    assert f2q(2.50000000000008, 0.001) == mpq(5,2)
    assert f2q(2.50000000000008, -50) == mpq(15637498706148,6254999482459)
    assert f2q(mpfr('0.500000011'), 1e-4) == mpq(1,2)
    assert f2q(mpfr('0.500000011'), 1e-5) == mpq(1,2)
    assert f2q(mpfr('0.500000011'), 1e-6) == mpq(1,2)
    assert f2q(mpfr('0.500000011'), 1e-7) == mpq(1,2)
    assert f2q(mpfr('0.500000011'), 1e-8) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-9) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-10) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-11) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-12) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-13) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-14) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-15) == mpq(22727273,45454545)
    assert f2q(mpfr('0.500000011'), 1e-16) == mpq(204545458,409090907)
    assert f2q(mpfr('0.500000011'), 1e-17) == mpq(204545458,409090907)


def test_get_emin_min():
    assert get_emin_min() in (-4611686018427387903, -1073741823)


def test_get_emax_max():
    assert get_emax_max() in (4611686018427387903, 1073741823)


def test_get_exp():
    assert get_exp(mpfr(5.232)) == 3

    pytest.raises(TypeError, lambda: get_exp(0))

    assert get_exp(mpfr('inf')) == 0
    assert get_exp(mpfr(0)) == 0


def test_set_exp():
    r = mpfr(4.55)

    assert set_exp(r, 4) == mpfr('9.0999999999999996')

    pytest.raises(TypeError, lambda: set_exp(r, mpz(4)))


def test_set_sign():
    r = mpfr(4.55)

    assert set_sign(r, False) == mpfr('4.5499999999999998')
    assert set_sign(r, True) == mpfr('-4.5499999999999998')

    pytest.raises(TypeError, lambda: set_sign(mpz(5), True))
    pytest.raises(TypeError, lambda: set_sign(r, 'oiio'))


def test_copy_sign():
    assert copy_sign(mpfr(4), mpfr(-2)) == mpfr('-4.0')

    pytest.raises(TypeError, lambda: copy_sign(mpfr(4), True))


def test_nan():
    x = nan()

    assert is_nan(x)


def test_inf():
    assert inf() == mpfr('inf')
    assert inf(-5) == mpfr('-inf')
    assert inf(mpz(-30)) == mpfr('-inf')

    pytest.raises(TypeError, lambda: inf(mpfr(30)))


def test_check_range():
    r = mpfr(4.55)

    assert check_range(r) == mpfr('4.5499999999999998')

    ctx = gmpy2.get_context()

    assert ctx.check_range(r) == mpfr('4.5499999999999998')

    pytest.raises(TypeError, lambda: ctx.check_range(mpz(5)))


def test_sign():
    a = mpq(3,11)

    assert sign(a) == 1
    assert sign(-a) == -1
    assert sign(mpq(0,5)) == 0

    pytest.raises(TypeError, lambda: sign('str'))

    a = mpfr("12.34")

    assert sign(-1.5) == -1
    assert sign(a) == 1
    assert sign(mpfr(0)) == 0
    assert sign(mpfr('inf')) == 1
    assert sign(mpfr('-inf')) == -1
    assert sign(mpfr('nan')) == 0
