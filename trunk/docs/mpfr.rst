Multiple-precision Reals
========================

gmpy2 replaces the *mpf* type from gmpy 1.x with a new *mpfr* type based on
the MPFR library. The new *mpfr* type supports correct rounding, selectable
rounding modes, and many trigonometric, exponential, and special functions. A
*context manager* is used to control precision, rounding modes, and the
behavior of exceptions.

The default precision of an *mpfr* is 53 bits - the same precision as Python's
*float* type. If the precison is changed, then ``mpfr(float('1.2'))`` differs
from ``mpfr('1.2')``. To take advantage of the higher precision provided by
the *mpfr* type, always pass constants as stringsl

::

    >>> import gmpy2
    >>> from gmpy2 import mpfr
    >>> mpfr('1.2')
    mpfr('1.2')
    >>> mpfr(float('1.2'))
    mpfr('1.2')
    >>> gmpy2.get_context().precision=100
    >>> mpfr('1.2')
    mpfr('1.2000000000000000000000000000006',100)
    >>> mpfr(float('1.2'))
    mpfr('1.1999999999999999555910790149937',100)
    >>>

Contexts
--------

A *context* is used to control the behavior of *mpfr* (and *mpc*) arithmetic.
In addition to controlling the precision, the rounding mode can be specified,
minimum and maximum exponent values can be changed, various exceptions can be
raised or ignored, gradual underflow can be enabled, and returning complex
results can be enabled.

``gmpy2.context()`` creates a new context with all options set to default.
``gmpy2.set_context(*ctx*)`` will set the active context to *ctx*.
``gmpy2.get_context()`` will return a reference to the active context. Note
that contexts are mutable: modifying the reference returned by get_context()
will modify the active context until a new context is enabled with
set_context().

The following example just modifies the precision. The remaining options will
be discussed later.

::

    >>> gmpy2.set_context(gmpy2.context())
    >>> gmpy2.get_context()
    context(precision=53, mpc_rprec=Default, mpc_iprec=Default,
            round=RoundToNearest, mpc_rround=Default, mpc_iround=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=False,
            trap_invalid=False, invalid=False,
            trap_erange=False, erange=False,
            trap_divzero=False, divzero=False,
            trap_expbound=False,
            allow_complex=False)
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
    >>>



