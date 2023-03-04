Multiple-precision Complex
==========================

.. currentmodule:: gmpy2

gmpy2 adds a multiple-precision complex type called `mpc` that is based
on the MPC library. The context manager settings for `mpfr` arithmetic are
applied to `mpc` arithmetic by default. It is possible to specify
different precision and rounding modes for both the real and imaginary
components of an `mpc`.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpc
    >>> gmpy2.set_context(gmpy2.context())
    >>> gmpy2.sqrt(mpc("1+2j"))
    mpc('1.272019649514069+0.78615137775742328j')
    >>> gmpy2.set_context(gmpy2.context(real_prec=60,imag_prec=70))
    >>> gmpy2.get_context()
    context(precision=53, real_prec=60, imag_prec=70,
            round=RoundToNearest, real_round=Default, imag_round=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=False,
            trap_invalid=False, invalid=False,
            trap_erange=False, erange=False,
            trap_divzero=False, divzero=False,
            allow_complex=False,
            rational_division=False,
            allow_release_gil=False)
    >>> gmpy2.sqrt(mpc("1+2j"))
    mpc('1.272019649514068965+0.78615137775742328606947j',(60,70))
    >>> gmpy2.set_context(gmpy2.context())

Exceptions are normally raised in Python when the result of a real operation
is not defined over the reals; for example, ``sqrt(-4)`` will raise an
exception. The default context in gmpy2 implements the same behavior but by
setting allow_complex to True, complex results will be returned.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpc
    >>> gmpy2.sqrt(-4)
    mpfr('nan')
    >>> gmpy2.get_context().allow_complex=True
    >>> gmpy2.sqrt(-4)
    mpc('0.0+2.0j')

The `mpc` type supports the `~mpc.__format__()` special method to
allow custom output formatting.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpc
    >>> a=gmpy2.sqrt(mpc("1+2j"))
    >>> a
    mpc('1.272019649514069+0.78615137775742328j')
    >>> "{0:.4.4Mf}".format(a)
    '(1.2720 0.7862)'
    >>> "{0:.4.4f}".format(a)
    '1.2720+0.7862j'
    >>> "{0:^20.4.4U}".format(a)
    '   1.2721+0.7862j   '
    >>> "{0:^20.4.4D}".format(a)
    '   1.2720+0.7861j   '

mpc Type
--------

.. autoclass:: mpc
   :special-members: __format__

mpc Functions
-------------

.. autofunction:: acos
.. autofunction:: acosh
.. autofunction:: asin
.. autofunction:: asinh
.. autofunction:: atan
.. autofunction:: atanh
.. autofunction:: cos
.. autofunction:: cosh
.. autofunction:: div_2exp
.. autofunction:: exp
.. autofunction:: is_nan
.. autofunction:: is_zero
.. autofunction:: log
.. autofunction:: log10
.. autofunction:: mpc_random
.. autofunction:: mul_2exp
.. autofunction:: norm
.. autofunction:: phase
.. autofunction:: polar
.. autofunction:: proj
.. autofunction:: rect
.. autofunction:: root_of_unity
.. autofunction:: sin
.. autofunction:: sin_cos
.. autofunction:: sinh
.. autofunction:: sqrt
.. autofunction:: tan
.. autofunction:: tanh
