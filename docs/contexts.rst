Contexts
========

.. currentmodule:: gmpy2

A `context` type is used to control the behavior
of `mpfr` and `mpc` arithmetic.  In addition to controlling the
precision, the rounding mode can be specified, minimum and maximum exponent
values can be changed, various exceptions can be raised or ignored, gradual
underflow can be enabled, and returning complex results can be enabled.

`context()` creates a new context with all options set to default.
`set_context()` will set the active context.  `get_context()` will
return a reference to the active context. Note that contexts are mutable:
modifying the reference returned by `get_context()` will modify the active
context until a new context is enabled with `set_context()`. The
`~context.copy()` method of a context will return a copy of the context.

The following example just modifies the precision. The remaining options will
be discussed later.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import mpfr
    >>> gmpy2.set_context(gmpy2.context())
    >>> gmpy2.get_context()
    context(precision=53, real_prec=Default, imag_prec=Default,
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
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997898')
    >>> gmpy2.get_context().precision=100
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997896964091736687316',100)
    >>> gmpy2.get_context().precision+=20
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997896964091736687312762351',120)
    >>> ctx=gmpy2.get_context()
    >>> ctx.precision+=20
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997896964091736687312762354406182',140)
    >>> gmpy2.set_context(gmpy2.context())
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997898')
    >>> ctx.precision+=20
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997898')
    >>> gmpy2.set_context(ctx)
    >>> gmpy2.sqrt(5)
    mpfr('2.2360679774997896964091736687312762354406183596116',160)
    >>> gmpy2.set_context(gmpy2.context())

Context Type
------------

.. autoclass:: context

Context Functions
-----------------

.. autofunction:: get_context
.. autofunction:: ieee
.. autofunction:: local_context
.. autofunction:: set_context

Contexts and the with statement
-------------------------------

Contexts can also be used in conjunction with Python's :keyword:`with`
statement to temporarily change the context settings for a block of code and
then restore the original settings when the block of code exits.

`local_context` first save the current context and then creates a new
context based on a context passed as the first argument, or the current context
if no context is passed. The new context is modified if any optional keyword
arguments are given. The original active context is restored when the block
completes.

In the following example, the current context is saved by `local_context`
and then the block begins with a copy of the default context and the precision
set to 100. When the block is finished, the original context is restored.

.. doctest::

    >>> print(gmpy2.sqrt(2))
    1.4142135623730951
    >>> with gmpy2.local_context(gmpy2.context(), precision=100) as ctx:
    ...   print(gmpy2.sqrt(2))
    ...   ctx.precision += 100
    ...   print(gmpy2.sqrt(2))
    ...
    1.4142135623730950488016887242092
    1.4142135623730950488016887242096980785696718753769480731766796
    >>> print(gmpy2.sqrt(2))
    1.4142135623730951

Contexts that implement the standard *single*, *double*, and *quadruple*
precision floating point types can be created using `ieee`.
