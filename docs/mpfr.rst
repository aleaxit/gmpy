Multiple-precision Reals
========================

The *mpfr* type is based on the MPFR library. The new *mpfr* type supports
correct rounding, selectable rounding modes, and many trigonometric,
exponential, and special functions. A *context manager* is used to control
precision, rounding modes, and the behavior of exceptions.

The default precision of an *mpfr* is 53 bits - the same precision as Python's
*float* type. If the precision is changed, then ``mpfr(float('1.2'))`` differs
from ``mpfr('1.2')``. To take advantage of the higher precision provided by
the *mpfr* type, always pass constants as strings.

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

A *context* is used to control the behavior of *mpfr* and *mpc* arithmetic.
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

::

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
            trap_expbound=False,
            allow_complex=False,
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
    >>>

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

.. note::
    It is possible to change the values of emin/emax such that previous *mpfr*
    values are no longer valid numbers but should either underflow to +/-0.0 or
    overflow to +/-Infinity. To raise an exception if this occurs, see
    **trap_expbound**.

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

    ::

        >>> gmpy2.set_context(gmpy2.context())
        >>> mpfr('nan') == mpfr('nan')
        False
        >>> gmpy2.get_context().trap_erange=True
        >>> mpfr('nan') == mpfr('nan')
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        gmpy2.RangeError: comparison with NaN
        >>>

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

**trap_expbound**
    This attribute controls whether or not an ``ExponentOutOfBoundsError``
    exception is raised if exponents in an operand are outside the current
    emin/emax limits.

**allow_complex**
    This attribute controls whether or not an *mpc* result can be returned if
    an *mpfr* result would normally not be possible.

Context Methods
---------------

**clear_flags()**
    Clear the underflow, overflow, inexact, invalid, erange, and divzero flags.

**copy()**
    Return a copy of the context.

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

::

    >>> with gmpy2.local_context(gmpy2.context(), precision=100) as ctx:
    ...   print(gmpy2.sqrt(2))
    ...   ctx.precision += 100
    ...   print(gmpy2.sqrt(2))
    ...
    1.4142135623730950488016887242092
    1.4142135623730950488016887242096980785696718753769480731766796
    >>>

A context object can also be used directly to create a context manager block.
However, instead of restoring the context to the active context when the
``with ...`` statement is executed, the restored context is the context used
before any keyword argument modifications.

The code:

::
    with gmpy2.ieee(64) as ctx:

is equivalent to:

::
    gmpy2.set_context(gmpy2.ieee(64))
    with gmpy2.local_context() as ctx:

Contexts that implement the standard *single*, *double*, and *quadruple* precision
floating point types can be created using **ieee()**.


mpfr Methods
------------

**as_integer_ratio()**
    Returns a 2-tuple containing the numerator and denominator after converting
    the *mpfr* object into the exact rational equivalent. The return 2-tuple
    is equivalent to Python's as_integer_ratio() method of built-in float
    objects.

**as_mantissa_exp()**
    Returns a 2-tuple containing the mantissa and exponent.

**as_simple_fraction()**
    Returns an *mpq* containing the simplest rational value that approximates
    the *mpfr* value with an error less than 1/(2**precision).

**conjugate()**
    Returns the complex conjugate. For *mpfr* objects, returns a copy of the
    original object.

**digits()**
    Returns a 3-tuple containing the mantissa, the exponent, and the number
    of bits of precision. The mantissa is represented as a string in the
    specified base with up to 'prec' digits. If 'prec' is 0, as many digits
    that are available are returned. No more digits than available given x's
    precision are returned. 'base' must be between 2 and 62, inclusive.

**is_integer()**
    Returns True if the *mpfr* object is an integer.

mpfr Attributes
---------------

**imag**
    Returns the imaginary component. For *mpfr* objects, returns 0.

**precision**
    Returns the precision of the *mpfr* object.

**rc**
    The result code (also known as ternary value in the MPFR documentation)
    is 0 if the value of the *mpfr* object is exactly equal to the exact,
    infinite precision value. If the result code is 1, then the value of the
    *mpfr* object is greater than the exact value. If the result code is -1,
    then the value of the *mpfr* object is less than the exact, infinite
    precision value.

**real**
    Returns the real component. For *mpfr* objects, returns a copy of the
    original object.

mpfr Functions
--------------

**acos(...)**
    acos(x) returns the arc-cosine of x. x is measured in radians. If
    context.allow_complex is True, then an *mpc* result will be returned for
    abs(x) > 1.

**acosh(...)**
    acosh(x) returns the inverse hyperbolic cosine of x.

**add(...)**
    add(x, y) returns x + y. The type of the result is based on the types of
    the arguments.

**agm(...)**
    agm(x, y) returns the arithmetic-geometric mean of x and y.

**ai(...)**
    ai(x) returns the Airy function of x.

**asin(...)**
    asin(x) returns the arc-sine of x. x is measured in radians. If
    context.allow_complex is True, then an *mpc* result will be returned for
    abs(x) > 1.

**asinh(...)**
    asinh(x) return the inverse hyperbolic sine of x.

**atan(...)**
    atan(x) returns the arc-tangent of x. x is measured in radians.

**atan2(...)**
    atan2(y, x) returns the arc-tangent of (y/x).

**atanh(...)**
    atanh(x) returns the inverse hyperbolic tangent of x. If
    context.allow_complex is True, then an *mpc* result will be returned for
    abs(x) > 1.

**cbrt(...)**
    cbrt(x) returns the cube root of x.

**ceil(...)**
    ceil(x) returns the 'mpfr' that is the smallest integer >= x.

**check_range(...)**
    check_range(x) return a new 'mpfr' with exponent that lies within the
    current range of emin and emax.

**const_catalan(...)**
    const_catalan([precision=0]) returns the Catalan's constant using the
    specified precision. If no precision is specified, the default precision
    is used.

**const_euler(...)**
    const_euler([precision=0]) returns the Euler's constant using the specified
    precision. If no precision is specified, the default precision is used.

**const_log2(...)**
    const_log2([precision=0]) returns the log2 constant using the specified
    precision. If no precision is specified, the default precision is used.

**const_pi(...)**
    const_pi([precision=0]) returns the constant pi using the specified
    precision. If no precision is specified, the default precision is used.

**context(...)**
    context() returns a new context manager controlling MPFR and MPC
    arithmetic.

**cos(...)**
    cos(x) returns the cosine of x. x is measured in radians.

**cosh(...)**
    cosh(x) returns the hyperbolic cosine of x.

**cot(...)**
    cot(x) returns the cotangent of x. x is measured in radians.

**coth(...)**
    coth(x) returns the hyperbolic cotangent of x.

**csc(...)**
    csc(x) returns the cosecant of x. x is measured in radians.

**csch(...)**
    csch(x) returns the hyperbolic cosecant of x.

**degrees(...)**
    degrees(x) converts an angle measurement x from radians to degrees.

**digamma(...)**
    digamma(x) returns the digamma of x.

**div(...)**
    div(x, y) returns x / y. The type of the result is based on the types of
    the arguments.

**div_2exp(...)**
    div_2exp(x, n) returns an 'mpfr' or 'mpc' divided by 2**n.

**eint(...)**
    eint(x) returns the exponential integral of x.

**erf(...)**
    erf(x) returns the error function of x.

**erfc(...)**
    erfc(x) returns the complementary error function of x.

**exp(...)**
    exp(x) returns e**x.

**exp10(...)**
    exp10(x) returns 10**x.

**exp2(...)**
    exp2(x) returns 2**x.

**expm1(...)**
    expm1(x) returns e**x - 1. expm1() is more accurate than exp(x) - 1 when
    x is small.

**f2q(...)**
    f2q(x[,err]) returns the simplest *mpq* approximating x to within relative
    error err. Default is the precision of x. Uses Stern-Brocot tree to find
    the simplest approximation. An *mpz* is returned if the denominator
    is 1. If err<0, error sought is 2.0 ** err.

**factorial(...)**
    factorial(n) returns the floating-point approximation to the factorial
    of n.

    See fac(n) to get the exact integer result.

**floor(...)**
    floor(x) returns the 'mpfr' that is the smallest integer <= x.

**fma(...)**
    fma(x, y, z) returns correctly rounded result of (x * y) + z.

**fmma(...)**
    fmma(x, y, z, t) returns correctly rounded result of (x * y) + (z * t).
    Requires MPFR 4.

**fmms(...)**
    fmms(x, y, z, t) returns correctly rounded result of (x * y) - (z * t).
    Requires MPFR 4.

**fmod(...)**
    fmod(x, y) returns x - n*y where n is the integer quotient of x/y, rounded
    to 0.

**fms(...)**
    fms(x, y, z) returns correctly rounded result of (x * y) - z.

**frac(...)**
    frac(x) returns the fractional part of x.

**frexp(...)**
    frexp(x) returns a tuple containing the exponent and mantissa of x.

**fsum(...)**
    fsum(iterable) returns the accurate sum of the values in the iterable.

**gamma(...)**
    gamma(x) returns the gamma of x.

**get_exp(...)**
    get_exp(mpfr) returns the exponent of an *mpfr*. Returns 0 for NaN or
    Infinity and sets the erange flag and will raise an exception if trap_erange
    is set.

**hypot(...)**
    hypot(y, x) returns square root of (x**2 + y**2).

**ieee(...)**
    ieee(bitwidth) returns a context with settings for 32-bit (aka single),
    64-bit (aka double), or 128-bit (aka quadruple) precision floating
    point types.

**inf(...)**
    inf(n) returns an *mpfr* initialized to Infinity with the same sign as n.
    If n is not given, +Infinity is returned.

**is_finite(...)**
    is_finite(x) returns True if x is an actual number (i.e. not NaN or
    Infinity).

**is_inf(...)**
    is_inf(x) returns True if x is Infinity or -Infinity.

    .. note::
        **is_inf()** is deprecated; please use **if_infinite()**.

**is_infinite(...)**
    is_infinite(x) returns True if x Infinity or -Infinity.

**is_nan(...)**
    is_nan(x) returns True if x is NaN (Not-A-Number).

**is_number(...)**
    is_number(x) returns True if x is an actual number (i.e. not NaN or
    Infinity).

    .. note::
        **is_number()** is deprecated; please use **is_finite()**.

**is_regular(...)**
    is_regular(x) returns True if x is not zero, NaN, or Infinity.

**is_signed(...)**
    is_signed(x) returns True if the sign bit of x is set.

**is_unordered(...)**
    is_unordered(x,y) returns True if either x and/or y is NaN.

**is_zero(...)**
    is_zero(x) returns True if x is zero.

**j0(...)**
    j0(x) returns the Bessel function of the first kind of order 0 of x.

**j1(...)**
    j1(x) returns the Bessel function of the first kind of order 1 of x.

**jn(...)**
    jn(x,n) returns the Bessel function of the first kind of order n of x.

**lgamma(...)**
    lgamma(x) returns a tuple containing the logarithm of the absolute value of
    gamma(x) and the sign of gamma(x)

**li2(...)**
    li2(x) returns the real part of dilogarithm of x.

**lngamma(...)**
    lngamma(x) returns the logarithm of gamma(x).

**log(...)**
    log(x) returns the natural logarithm of x.

**log10(...)**
    log10(x) returns the base-10 logarithm of x.

**log1p(...)**
    log1p(x) returns the natural logarithm of (1+x).

**log2(...)**
    log2(x) returns the base-2 logarithm of x.

**max2(...)**
    max2(x, y) returns the maximum of x and y. The result may be rounded to
    match the current context. Use the builtin max() to get an exact copy of
    the largest object without any rounding.

**min2(...)**
    min2(x, y) returns the minimum of x and y. The result may be rounded to
    match the current context. Use the builtin min() to get an exact copy of
    the smallest object without any rounding.

**modf(...)**
    modf(x) returns a tuple containing the integer and fractional portions
    of x.

**mpfr(...)**
    mpfr() returns and *mpfr* object set to 0.0.

    mpfr(n[, precision=0]) returns an *mpfr* object after converting a numeric
    value n. If no precision, or a precision of 0, is specified; the precision
    is taken from the current context.

    mpfr(s[, precision=0[, [base=0]]) returns an *mpfr* object after converting
    a string 's' made up of digits in the given base, possibly with fractional
    part (with period as a separator) and/or exponent (with exponent marker
    'e' for base<=10, else '@'). If no precision, or a precision of 0, is
    specified; the precision is taken from the current context. The base of the
    string representation must be 0 or in the interval 2 ... 62. If the base
    is 0, the leading digits of the string are used to identify the base: 0b
    implies base=2, 0x implies base=16, otherwise base=10 is assumed.

**mpfr_from_old_binary(...)**
    mpfr_from_old_binary(string) returns an *mpfr* from a GMPY 1.x binary mpf
    format. Please use to_binary()/from_binary() to convert GMPY2 objects to or
    from a binary format.

**mpfr_grandom(...)**
    mpfr_grandom(random_state) returns two random numbers with Gaussian
    distribution. The parameter *random_state* must be created by random_state()
    first.

**mpfr_random(...)**
    mpfr_random(random_state) returns a uniformly distributed number between
    [0,1]. The parameter *random_state* must be created by random_state() first.

**mul(...)**
    mul(x, y) returns x * y. The type of the result is based on the types of
    the arguments.

**mul_2exp(...)**
    mul_2exp(x, n) returns 'mpfr' or 'mpc' multiplied by 2**n.

**nan(...)**
    nan() returns an 'mpfr' initialized to NaN (Not-A-Number).

**next_above(...)**
    next_above(x) returns the next 'mpfr' from x toward +Infinity.

**next_below(...)**
    next_below(x) returns the next 'mpfr' from x toward -Infinity.

**radians(...)**
    radians(x) converts an angle measurement x from degrees to radians.

**rec_sqrt(...)**
    rec_sqrt(x) returns the reciprocal of the square root of x.

**reldiff(...)**
    reldiff(x, y) returns the relative difference between x and y. Result is
    equal to abs(x-y)/x.

**remainder(...)**
    remainder(x, y) returns x - n*y where n is the integer quotient of x/y,
    rounded to the nearest integer and ties rounded to even.

**remquo(...)**
    remquo(x, y) returns a tuple containing the remainder(x,y) and the low bits
    of the quotient.

**rint(...)**
    rint(x) returns x rounded to the nearest integer using the current rounding
    mode.

**rint_ceil(...)**
    rint_ceil(x) returns x rounded to the nearest integer by first rounding to
    the next higher or equal integer and then, if needed, using the current
    rounding mode.

**rint_floor(...)**
    rint_floor(x) returns x rounded to the nearest integer by first rounding to
    the next lower or equal integer and then, if needed, using the current
    rounding mode.

**rint_round(...)**
    rint_round(x) returns x rounded to the nearest integer by first rounding to
    the nearest integer (ties away from 0) and then, if needed, using the
    current rounding mode.

**rint_trunc(...)**
    rint_trunc(x) returns x rounded to the nearest integer by first rounding
    towards zero and then, if needed, using the current rounding mode.

**root(...)**
    root(x, n) returns n-th root of x. The result always an *mpfr*.

**round2(...)**
    round2(x[, n]) returns x rounded to n bits. Uses default precision if n is
    not specified. See round_away() to access the mpfr_round() function. Use
    the builtin round() to round x to n decimal digits.

**round_away(...)**
    round_away(x) returns an *mpfr* by rounding x the nearest integer, with
    ties rounded away from 0.

**sec(...)**
    sec(x) returns the secant of x. x is measured in radians.

**sech(...)**
    sech(x) returns the hyperbolic secant of x.

**set_exp(...)**
    set_exp(x, n) sets the exponent of a given *mpfr* to n. If n is outside the
    range of valid exponents, set_exp() will set the erange flag and either
    return the original value or raise an exception if trap_erange is set.

**set_sign(...)**
    set_sign(x, bool) returns a copy of x with it's sign bit set if *bool*
    evaluates to True.

**sign(...)**
    sign(x) returns -1 if x < 0, 0 if x == 0, or +1 if x >0.

**sin(...)**
    sin(x) returns the sine of x. x is measured in radians.

**sin_cos(...)**
    sin_cos(x) returns a tuple containing the sine and cosine of x. x is
    measured in radians.

**sinh(...)**
    sinh(x) returns the hyberbolic sine of x.

**sinh_cosh(...)**
    sinh_cosh(x) returns a tuple containing the hyperbolic sine and cosine of
    x.

**sqrt(...)**
    sqrt(x) returns the square root of x. If x is integer, rational, or real,
    then an *mpfr* will be returned. If x is complex, then an *mpc* will
    be returned. If context.allow_complex is True, negative values of x
    will return an *mpc*.

**square(...)**
    square(x) returns x * x. The type of the result is based on the types of
    the arguments.

**sub(...)**
    sub(x, y) returns x - y. The type of the result is based on the types of
    the arguments.

**tan(...)**
    tan(x) returns the tangent of x. x is measured in radians.

**tanh(...)**
    tanh(x) returns the hyperbolic tangent of x.

**trunc(...)**
    trunc(x) returns an 'mpfr' that is x truncated towards 0. Same as
    x.floor() if x>=0 or x.ceil() if x<0.

**y0(...)**
    y0(x) returns the Bessel function of the second kind of order 0 of x.

**y1(...)**
    y1(x) returns the Bessel function of the second kind of order 1 of x.

**yn(...)**
    yn(x,n) returns the Bessel function of the second kind of order n of x.

**zero(...)**
    zero(n) returns an *mpfr* initialized to 0.0 with the same sign as n.
    If n is not given, +0.0 is returned.

**zeta(...)**
    zeta(x) returns the Riemann zeta of x.

mpfr Formatting
---------------

The *mpfr* type supports the __format__() special method to allow custom output
formatting.

**__format__(...)**
    x.__format__(fmt) returns a Python string by formatting 'x' using the
    format string 'fmt'. A valid format string consists of:

    |     optional alignment code:
    |        '<' -> left shifted in field
    |        '>' -> right shifted in field
    |        '^' -> centered in field
    |     optional leading sign code
    |        '+' -> always display leading sign
    |        '-' -> only display minus for negative values
    |        ' ' -> minus for negative values, space for positive values
    |     optional width.precision
    |     optional rounding mode:
    |        'U' -> round toward plus infinity
    |        'D' -> round toward minus infinity
    |        'Y' -> round away from zero
    |        'Z' -> round toward zero
    |        'N' -> round to nearest
    |     optional conversion code:
    |        'a','A' -> hex format
    |        'b'     -> binary format
    |        'e','E' -> scientific format
    |        'f','F' -> fixed point format
    |        'g','G' -> fixed or scientific format

    .. note::
        The formatting codes must be specified in the order shown above.

::

    >>> import gmpy2
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



