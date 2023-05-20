Multiple-precision Reals
========================

.. currentmodule:: gmpy2

The `mpfr` type is based on the MPFR library. The new `mpfr` type supports
correct rounding, selectable rounding modes, and many trigonometric,
exponential, and special functions. A *context manager* is used to control
precision, rounding modes, and the behavior of exceptions.

The default precision of an `mpfr` is 53 bits - the same precision as Python's
`float` type. If the precision is changed, then ``mpfr(float('1.2'))`` differs
from ``mpfr('1.2')``. To take advantage of the higher precision provided by
the `mpfr` type, always pass constants as strings.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpfr
    >>> gmpy2.set_context(gmpy2.context())
    >>> mpfr('1.2')
    mpfr('1.2')
    >>> mpfr(float('1.2'))
    mpfr('1.2')
    >>> gmpy2.get_context().precision=100
    >>> mpfr('1.2')
    mpfr('1.2000000000000000000000000000006',100)
    >>> mpfr(float('1.2'))
    mpfr('1.1999999999999999555910790149937',100)

The `mpfr` type supports the `~mpfr.__format__()` special method to
allow custom output formatting.

.. doctest::

    >>> from gmpy2 import mpfr
    >>> a=mpfr("1.23456")
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

mpfr Type
---------

.. autoclass:: mpfr
   :special-members: __format__

mpfr Functions
--------------

.. autofunction:: agm
.. autofunction:: ai
.. autofunction:: atan2
.. autofunction:: cbrt
.. autofunction:: ceil
.. autofunction:: check_range
.. autofunction:: cmp
.. autofunction:: const_catalan
.. autofunction:: const_euler
.. autofunction:: const_log2
.. autofunction:: const_pi
.. autofunction:: cot
.. autofunction:: coth
.. autofunction:: csc
.. autofunction:: csch
.. autofunction:: degrees
.. autofunction:: digamma
.. autofunction:: eint
.. autofunction:: erf
.. autofunction:: erfc
.. autofunction:: exp10
.. autofunction:: exp2
.. autofunction:: expm1
.. autofunction:: factorial
.. autofunction:: floor
.. autofunction:: fmma
.. autofunction:: fmms
.. autofunction:: fmod
.. autofunction:: frac
.. autofunction:: frexp
.. autofunction:: fsum
.. autofunction:: gamma
.. autofunction:: gamma_inc
.. autofunction:: get_exp
.. autofunction:: hypot
.. autofunction:: inf
.. autofunction:: is_finite
.. autofunction:: is_infinite
.. autofunction:: is_regular
.. autofunction:: is_signed
.. autofunction:: is_unordered
.. autofunction:: j0
.. autofunction:: j1
.. autofunction:: jn
.. autofunction:: lgamma
.. autofunction:: li2
.. autofunction:: lngamma
.. autofunction:: log1p
.. autofunction:: log2
.. autofunction:: maxnum
.. autofunction:: minnum
.. autofunction:: modf
.. autofunction:: mpfr_from_old_binary
.. autofunction:: mpfr_grandom
.. autofunction:: mpfr_nrandom
.. autofunction:: mpfr_random
.. autofunction:: nan
.. autofunction:: next_above
.. autofunction:: next_below
.. autofunction:: radians
.. autofunction:: rec_sqrt
.. autofunction:: reldiff
.. autofunction:: remainder
.. autofunction:: remquo
.. autofunction:: rint
.. autofunction:: rint_ceil
.. autofunction:: rint_floor
.. autofunction:: rint_round
.. autofunction:: rint_trunc
.. autofunction:: root
.. autofunction:: rootn
.. autofunction:: round2
.. autofunction:: round_away
.. autofunction:: sec
.. autofunction:: sech
.. autofunction:: set_exp
.. autofunction:: set_sign
.. autofunction:: sign
.. autofunction:: sinh_cosh
.. autofunction:: trunc
.. autofunction:: y0
.. autofunction:: y1
.. autofunction:: yn
.. autofunction:: zero
.. autofunction:: zeta
.. autofunction:: get_max_precision
.. autofunction:: get_emax_max
.. autofunction:: get_emin_min
.. autofunction:: copy_sign
.. autofunction:: can_round
.. autofunction:: free_cache
