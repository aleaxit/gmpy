Overview
========

.. currentmodule:: gmpy2

Tutorial
--------

The `mpz` type is compatible with Python's built-in `int` type but is
significantly faster for large values. The cutover point for performance
varies, but can be as low as 20 to 40 digits. A variety of additional integer
functions are provided.

Operator overloading is fully supported. Coversion from native Python types is
optimized for performance.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpz,mpq,mpfr,mpc
    >>> gmpy2.set_context(gmpy2.context())
    >>> mpz(99) * 43
    mpz(4257)
    >>> pow(mpz(99), 37, 59)
    mpz(18)
    >>> gmpy2.isqrt(99)
    mpz(9)
    >>> gmpy2.isqrt_rem(99)
    (mpz(9), mpz(18))
    >>> gmpy2.gcd(123,27)
    mpz(3)
    >>> gmpy2.lcm(123,27)
    mpz(1107)
    >>> (mpz(123) + 12) / 5
    mpfr('27.0')
    >>> (mpz(123) + 12) // 5
    mpz(27)
    >>> (mpz(123) + 12) / 5.0
    mpfr('27.0')

The `mpq` type is compatible with the `~fractions.Fraction` type included with
Python.

.. doctest::

    >>> mpq(3,7)/7
    mpq(3,49)
    >>> mpq(45,3) * mpq(11,8)
    mpq(165,8)

gmpy2 supports correctly rounded arbitrary precision real and complex arithmetic
via the MPFR and MPC libraries. Floating point contexts are used to control precision,
rounding modes, and exceptional conditions. For example, division by zero can either
return an Infinity or raise an exception.

.. doctest::

    >>> gmpy2.set_context(gmpy2.context())
    >>> mpfr(1)/7
    mpfr('0.14285714285714285')
    >>> gmpy2.get_context().precision=200
    >>> mpfr(1)/7
    mpfr('0.1428571428571428571428571428571428571428571428571428571428571',200)
    >>> gmpy2.get_context()
    context(precision=200, real_prec=Default, imag_prec=Default,
            round=RoundToNearest, real_round=Default, imag_round=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=True,
            trap_invalid=False, invalid=False,
            trap_erange=False, erange=False,
            trap_divzero=False, divzero=False,
            allow_complex=False,
            rational_division=False,
            allow_release_gil=False)
    >>> mpfr(1)/0
    mpfr('inf')
    >>> gmpy2.get_context().trap_divzero=True
    >>> mpfr(1)/0  # doctest: +IGNORE_EXCEPTION_DETAIL +ELLIPSIS
    Traceback (most recent call last):
    ...
    gmpy2.DivisionByZeroError: division by zero
    >>> gmpy2.get_context()
    context(precision=200, real_prec=Default, imag_prec=Default,
            round=RoundToNearest, real_round=Default, imag_round=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=True,
            trap_invalid=False, invalid=False,
            trap_erange=False, erange=False,
            trap_divzero=True, divzero=True,
            allow_complex=False,
            rational_division=False,
            allow_release_gil=False)
    >>> gmpy2.sqrt(mpfr(-2))
    mpfr('nan')
    >>> gmpy2.get_context().allow_complex=True
    >>> gmpy2.get_context().precision=53
    >>> gmpy2.sqrt(mpfr(-2))
    mpc('0.0+1.4142135623730951j')
    >>>
    >>> gmpy2.set_context(gmpy2.context())
    >>> with gmpy2.local_context() as ctx:
    ...   print(gmpy2.const_pi())
    ...   ctx.precision+=20
    ...   print(gmpy2.const_pi())
    ...   ctx.precision+=20
    ...   print(gmpy2.const_pi())
    ...
    3.1415926535897931
    3.1415926535897932384628
    3.1415926535897932384626433831
    >>> print(gmpy2.const_pi())
    3.1415926535897931


Miscellaneous gmpy2 Functions
-----------------------------

.. autofunction:: digits
.. autofunction:: from_binary
.. autofunction:: license
.. autofunction:: mp_limbsize
.. autofunction:: mp_version
.. autofunction:: mpc_version
.. autofunction:: mpfr_version
.. autofunction:: random_state
.. autofunction:: to_binary
.. autofunction:: version

Generic gmpy2 Functions
-----------------------

.. autofunction:: add
.. autofunction:: div
.. autofunction:: mul
.. autofunction:: sub

.. autofunction:: square

.. autofunction:: f2q

.. autofunction:: fma
.. autofunction:: fms

.. autofunction:: cmp_abs

Exceptions
----------

.. autoexception:: RangeError
.. autoexception:: InexactResultError
.. autoexception:: OverflowResultError
.. autoexception:: UnderflowResultError
.. autoexception:: InvalidOperationError
.. autoexception:: DivisionByZeroError
