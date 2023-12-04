Tutorial
========

.. currentmodule:: gmpy2

Start by importing the contents of the package with:

.. doctest::

    >>> from gmpy2 import *

.. note::

    The use of ``from gmpy2 import *`` is not recommended in real code.  The
    names in gmpy2 have been chosen to avoid conflict with Python's builtin
    names but gmpy2 does use names that may conflict with other modules or
    variable names.  In normal usage you’ll probably only want to import the
    classes and functions that you actually need.

Lets look first on some examples of arbitrary precision arithmetic with
integer and rational types:

.. doctest::

    >>> mpz(99) * 43
    mpz(4257)
    >>> pow(mpz(99), 37, 59)
    mpz(18)
    >>> isqrt(99)
    mpz(9)
    >>> isqrt_rem(99)
    (mpz(9), mpz(18))
    >>> gcd(123, 27)
    mpz(3)
    >>> lcm(123, 27)
    mpz(1107)
    >>> (mpz(123) + 12) / 5
    mpfr('27.0')
    >>> (mpz(123) + 12) // 5
    mpz(27)
    >>> (mpz(123) + 12) / 5.0
    mpfr('27.0')
    >>> mpz('123') + 1
    mpz(124)
    >>> 10 - mpz(1)
    mpz(9)
    >>> is_prime(17)
    True
    >>> mpz('1_000_000')
    mpz(1000000)
    >>> mpq(3, 7)/7
    mpq(3,49)
    >>> mpq(45, 3) * mpq(11, 8)
    mpq(165,8)
    >>> mpq(1, 7) * 11
    mpq(11,7)

But gmpy2 also supports correctly rounded multiple precision real and complex
arithmetic.  The following example shows how to control precision settings and
rounding modes:

.. doctest::

    >>> mpfr('1.2')
    mpfr('1.2')
    >>> mpfr(float('1.2'))
    mpfr('1.2')
    >>> ctx = get_context()
    >>> ctx.precision
    53
    >>> ctx.precision = 100
    >>> mpfr('1.2')
    mpfr('1.2000000000000000000000000000006',100)
    >>> mpfr(float('1.2'))
    mpfr('1.1999999999999999555910790149937',100)
    >>> ctx.precision = 53
    >>> ctx.round = RoundUp
    >>> const_pi()
    mpfr('3.1415926535897936')
    >>> ctx.round = RoundToNearest
    >>> const_pi()
    mpfr('3.1415926535897931')

You have seen, that if the precision is changed, then ``mpfr(float('1.2'))``
differs from ``mpfr('1.2')``.  To take advantage of the higher precision
provided by the `mpfr` type, always pass constants as strings.

Floating point contexts also are used to control exceptional conditions.  For
example, division by zero can either return a floating-point positive infinity
(default) or raise an exception.

.. doctest::

    >>> ctx.divzero
    False
    >>> mpfr(1)/0
    mpfr('inf')
    >>> ctx.trap_divzero = True
    >>> mpfr(1)/0
    Traceback (most recent call last):
    ...
    gmpy2.DivisionByZeroError: division by zero
    >>> ctx.divzero
    True

Exceptions are normally raised in Python when the result of a real operation is
not defined over the reals; for example, ``math.sqrt(-2)`` will raise a
:exc:`ValueError` exception.  The default context in gmpy2 implements similar
behavior, but by setting :attr:`~context.allow_complex` flag, complex results
will be returned.

.. doctest::

    >>> sqrt(mpfr(-2))
    mpfr('nan')
    >>> ctx.allow_complex = True
    >>> sqrt(mpfr(-2))
    mpc('0.0+1.4142135623730951j')

Contexts can also be used as context managers in conjunction with Python's
:keyword:`with` statement to temporarily change the current context settings
for a block of code.

.. doctest::

    >>> print(const_pi())
    3.1415926535897931
    >>> with context(precision=100) as ctx:
    ...   print(const_pi())
    ...   ctx.precision += 20
    ...   print(const_pi())
    ...
    3.1415926535897932384626433832793
    3.1415926535897932384626433832795028847
    >>> print(const_pi())
    3.1415926535897931

It's possible to set different precision settings for real and imaginary
components.

.. doctest::

    >>> ctx = get_context()
    >>> ctx.real_prec = 60
    >>> ctx.imag_prec = 70
    >>> sqrt(mpc('1+2j'))
    mpc('1.272019649514068965+0.78615137775742328606947j',(60,70))

All gmpy2 numeric types support Python's "new style" string formatting
available in `formatted string literals
<https://docs.python.org/3/tutorial/inputoutput.html#tut-f-strings>`_ or with
:meth:`str.format`; see `Format Specification Mini-Language
<https://docs.python.org/3/library/string.html#formatspec>`_ for a description
of the standard formatting syntax.  The precision value optionally can be
followed by the rounding mode type ('U' to round toward plus infinity, 'D' to
round toward minus infinity, 'Y' to round away from zero, 'Z' to round toward
zero and 'N' - round to the nearest value.

.. doctest::

    >>> a = mpfr("1.23456")
    >>> "{0:15.3f}".format(a)
    '          1.235'
    >>> "{0:15.3Uf}".format(a)
    '          1.235'
    >>> "{0:15.3Df}".format(a)
    '          1.234'
    >>> "{0:.3Df}".format(a)
    '1.234'
    >>> "{0:+.3Df}".format(a)
    '+1.234'
