Contexts
========

.. currentmodule:: gmpy2

A *context* is used to control the behavior of :class:`mpfr` and :class:`mpc` arithmetic.
In addition to controlling the precision, the rounding mode can be specified,
minimum and maximum exponent values can be changed, various exceptions can be
raised or ignored, gradual underflow can be enabled, and returning complex
results can be enabled.

``gmpy2.context()`` creates a new context with all options set to default.
``gmpy2.set_context(ctx)`` will set the active context to *ctx*.
``gmpy2.get_context()`` will return a reference to the active context. Note
that contexts are mutable: modifying the reference returned by get_context()
will modify the active context until a new context is enabled with
set_context(). The ``copy()`` method of a context will return a copy of the
context.

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

.. autoclass:: _context
   :members:

Context Attributes
------------------

**precision**
    This attribute controls the precision of an *mpfr* result. The precision
    is specified in bits, not decimal digits. The maximum precision that can
    be specified is platform dependent and can be retrieved with
    **get_max_precision()**.

.. note::
    Specifying a value for precision that is too close to the maximum precision
    will cause the MPFR library to fail.

**real_prec**
    This attribute controls the precision of the real part of an *mpc* result.
    If the value is ``Default``, then the value of the precision attribute is
    used.

**imag_prec**
    This attribute controls the precision of the imaginary part of an *mpc*
    result. If the value is ``Default``, then the value of real_prec is used.

**round**
    There are five rounding modes available to *mpfr* types:

    ``RoundAwayZero``
        The result is rounded away from 0.0.

    ``RoundDown``
        The result is rounded towards -Infinity.

    ``RoundToNearest``
        Round to the nearest value; ties are rounded to an even value.

    ``RoundToZero``
        The result is rounded towards 0.0.

    ``RoundUp``
        The result is rounded towards +Infinity.

**real_round**
    This attribute controls the rounding mode for the real part of an *mpc*
    result. If the value is ``Default``, then the value of the round attribute
    is used. Note: ``RoundAwayZero`` is not a valid rounding mode for *mpc*.

**imag_round**
    This attribute controls the rounding mode for the imaginary part of an
    *mpc* result. If the value is ``Default``, then the value of the real_round
    attribute is used. Note: ``RoundAwayZero`` is not a valid rounding mode for
    *mpc*.

**emax**
    This attribute controls the maximum allowed exponent of an *mpfr* result.
    The maximum exponent is platform dependent and can be retrieved with
    **get_emax_max()**.

**emin**
    This attribute controls the minimum allowed exponent of an *mpfr* result.
    The minimum exponent is platform dependent and can be retrieved with
    **get_emin_min()**.

**subnormalize**
    The usual IEEE-754 floating point representation supports gradual underflow
    when the minimum exponent is reached. The MFPR library does not enable
    gradual underflow by default but it can be enabled to precisely mimic the
    results of IEEE-754 floating point operations.

**trap_underflow**
    If set to ``False``, a result that is smaller than the smallest possible
    *mpfr* given the current exponent range will be replaced by +/-0.0. If set
    to ``True``, an ``UnderflowResultError`` exception is raised.

**underflow**
    This flag is not user controllable. It is automatically set if a result
    underflowed to +/-0.0 and trap_underflow is ``False``.

**trap_overflow**
    If set to ``False``, a result that is larger than the largest possible
    *mpfr* given the current exponent range will be replaced by +/-Infinity. If
    set to ``True``, an ``OverflowResultError`` exception is raised.

**overflow**
    This flag is not user controllable. It is automatically set if a result
    overflowed to +/-Infinity and trap_overflow is ``False``.

**trap_inexact**
    This attribute controls whether or not an ``InexactResultError`` exception
    is raised if an inexact result is returned. To check if the result is
    greater or less than the exact result, check the **rc** attribute of the
    *mpfr* result.

**inexact**
    This flag is not user controllable. It is automatically set if an inexact
    result is returned.

**trap_invalid**
    This attribute controls whether or not an ``InvalidOperationError``
    exception is raised if a numerical result is not defined. A special
    NaN (Not-A-Number) value will be returned if an exception is not raised.
    The ``InvalidOperationError`` is a sub-class of Python's ``ValueError``.

    For example, ``gmpy2.sqrt(-2)`` will normally return *mpfr('nan')*.
    However, if allow_complex is set to ``True``, then an *mpc* result will
    be returned.

**invalid**
    This flag is not user controllable. It is automatically set if an invalid
    (Not-A-Number) result is returned.

**trap_erange**
    This attribute controls whether or not a ``RangeError`` exception is raised
    when certain operations are performed on NaN and/or Infinity values.
    Setting trap_erange to ``True`` can be used to raise an exception if
    comparisons are attempted with a NaN.

    .. doctest::

        >>> gmpy2.set_context(gmpy2.context())
        >>> mpfr('nan') == mpfr('nan')
        False
        >>> gmpy2.get_context().trap_erange=True
        >>> mpfr('nan') == mpfr('nan')
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        gmpy2.RangeError: comparison with NaN
        >>> gmpy2.set_context(gmpy2.context())

**erange**
    This flag is not user controllable. It is automatically set if an erange
    error occurred.

**trap_divzero**
    This attribute controls whether or not a ``DivisionByZeroError`` exception
    is raised if division by 0 occurs. The ``DivisionByZeroError`` is a
    sub-class of Python's ``ZeroDivisionError``.

**divzero**
    This flag is not user controllable. It is automatically set if a division
    by zero occurred and NaN result was returned.

**allow_complex**
    This attribute controls whether or not an *mpc* result can be returned if
    an *mpfr* result would normally not be possible.

**rational_division**
    If set to ``True``, *mpz* / *mpz* will return an *mpq* instead of an *mpfr*.

**allow_release_gil**
    If set to ``True``, many *mpz* and *mpq* computations will release the GIL.

    This is considered an experimental feature.

Context Methods
---------------

**abs**

**acos**

**acosh**

**add**

**agm**

**ai**

**asin**

**asinh**

**atan**

**atan2**

**atanh**

**cbrt**

**ceil**

**check_range**

**clear_flags()**
    Clear the underflow, overflow, inexact, invalid, erange, and divzero flags.

**const_catalan**

**const_euler**

**const_log**

**const_pi**

**copy()**
    Return a copy of the context.

**cos**

**cosh**

**cot**

**coth**

**csc**

**degrees**

**digamma**

**div**

**div_2exp**

**divmod**

**eint**

**erf**

**erfc**

**exp**

**exp10**

**exp2**

**expm1**

**factorial**

**floor**

**floor_div**

**fma**

**fmma**

**fmms**

**fmod**

**fms**

**frac**

**frexp**

**fsum**

**gamma**

**hypot**

**is_finite**

**is_infinite**

**is_integer**

**is_nan**

**is _regular**

**is_signed**

**is_zero**

**j0**

**j1**

**jn**

**lgamma**

**li2**

**lngamma**

**log**

**log10**

**log1p**

**log2**

**maxnum**

**minnum**

**minus**

**mod**

**modf**

**mul**

**mul_2exp**

**next_above**

**next_below**

**next_toward**

**norm**

**phase**

**plus**

**polar**

**pow**

**proj**

**radians**

**rec_sqrt**

**rect**

**reldiff**

**remainder**

**remquo**

**rint**

**rint_ceil**

**rint_floor**

**rint_round**

**rint_trunc**

**root**

**root_of_unity**

**rootn**

**round**

**round2**

**round_away**

**sec**

**sech**

**sin**

**sin_cos**

**sinh**

**sinh_cosh**

**sqrt**

**square**

**sub**

**subnormalize**

**tan**

**tanh**

**trunc**

**y0**

**y1**

**yn**

**zeta**

Contexts and the with statement
-------------------------------

Contexts can also be used in conjunction with Python's ``with ...`` statement to
temporarily change the context settings for a block of code and then restore the
original settings when the block of code exits.

``gmpy2.local_context()`` first save the current context and then creates a new
context based on a context passed as the first argument, or the current context
if no context is passed. The new context is modified if any optional keyword
arguments are given. The original active context is restored when the block
completes.

In the following example, the current context is saved by ``gmpy2.local_context()``
and then the block begins with a copy of the default context and the precision
set to 100. When the block is finished, the original context is restored.

.. doctest::

    >>> with gmpy2.local_context(gmpy2.context(), precision=100) as ctx:
    ...   print(gmpy2.sqrt(2))
    ...   ctx.precision += 100
    ...   print(gmpy2.sqrt(2))
    ...
    1.4142135623730950488016887242092
    1.4142135623730950488016887242096980785696718753769480731766796

A context object can also be used directly to create a context manager block.
However, instead of restoring the context to the active context when the
``with ...`` statement is executed, the restored context is the context used
before any keyword argument modifications.

The code:

.. code-block:: python

    with gmpy2.ieee(64) as ctx:

is equivalent to:

.. code-block:: python

    gmpy2.set_context(gmpy2.ieee(64))
    with gmpy2.local_context() as ctx:

Contexts that implement the standard *single*, *double*, and *quadruple* precision
floating point types can be created using **ieee()**.
