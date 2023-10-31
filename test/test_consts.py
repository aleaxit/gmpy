import pytest

import gmpy2
from gmpy2 import const_catalan, const_euler, const_log2, const_pi, mpfr


def test_pi():
    assert const_pi() == mpfr('3.1415926535897931')
    assert const_pi(100) == mpfr('3.1415926535897932384626433832793',100)

    pytest.raises(TypeError, lambda: const_pi(100,200))
    pytest.raises(TypeError, lambda: const_pi(prec=100))

    assert const_pi(precision=100) == mpfr('3.1415926535897932384626433832793',100)
    assert gmpy2.ieee(32).const_pi() == mpfr('3.14159274',24)

    pytest.raises(TypeError, lambda: gmpy2.ieee(32).const_pi(100))

    assert gmpy2.ieee(128).const_pi() == mpfr('3.1415926535897932384626433832795028',113)


def test_log2():
    assert const_log2() == mpfr('0.69314718055994529')
    assert const_log2(100) == mpfr('0.69314718055994530941723212145798',100)

    pytest.raises(TypeError, lambda: const_log2(100,200))
    pytest.raises(TypeError, lambda: const_log2(prec=100))

    assert const_log2(precision=100) == mpfr('0.69314718055994530941723212145798',100)
    assert gmpy2.ieee(32).const_log2() == mpfr('0.693147182',24)

    pytest.raises(TypeError, lambda: gmpy2.ieee(32).const_log2(100))

    assert gmpy2.ieee(128).const_log2() == mpfr('0.693147180559945309417232121458176575',113)


def test_catalan():
    assert const_catalan() == mpfr('0.91596559417721901')
    assert const_catalan(100) == mpfr('0.91596559417721901505460351493252',100)

    pytest.raises(TypeError, lambda: const_catalan(100,200))
    pytest.raises(TypeError, lambda: const_catalan(prec=100))

    assert gmpy2.const_catalan(precision=100) == mpfr('0.91596559417721901505460351493252',100)
    assert gmpy2.ieee(32).const_catalan() == mpfr('0.915965617',24)

    pytest.raises(TypeError, lambda: gmpy2.ieee(32).const_catalan(100))

    assert gmpy2.ieee(128).const_catalan() == mpfr('0.915965594177219015054603514932384146',113)


def test_euler():
    assert const_euler() == mpfr('0.57721566490153287')
    assert const_euler(100) == mpfr('0.57721566490153286060651209008234',100)

    pytest.raises(TypeError, lambda: const_euler(100,200))
    pytest.raises(TypeError, lambda: const_euler(prec=100))

    assert const_euler(precision=100) == mpfr('0.57721566490153286060651209008234',100)
    assert gmpy2.ieee(32).const_euler() == mpfr('0.577215672',24)

    pytest.raises(TypeError, lambda: gmpy2.ieee(32).const_euler(100))

    assert gmpy2.ieee(128).const_euler() == mpfr('0.577215664901532860606512090082402471',113)
