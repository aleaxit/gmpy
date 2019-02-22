Multiple-precision Complex
==========================

gmpy2 adds a multiple-precision complex type called *mpc* that is based on the
MPC library. The context manager settings for *mpfr* arithmetic are applied to
*mpc* arithmetic by default. It is possible to specify different precision and
rounding modes for both the real and imaginary components of an *mpc*.

::

    >>> import gmpy2
    >>> from gmpy2 import mpc
    >>> gmpy2.sqrt(mpc("1+2j"))
    mpc('1.272019649514069+0.78615137775742328j')
    >>> gmpy2.get_context(real_prec=100,imag_prec=200)
    context(precision=53, real_prec=100, imag_prec=200,
            round=RoundToNearest, real_round=Default, imag_round=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=True,
            trap_invalid=False, invalid=False,
            trap_erange=False, erange=False,
            trap_divzero=False, divzero=False,
            trap_expbound=False,
            allow_complex=False)
    >>> gmpy2.sqrt(mpc("1+2j"))
    mpc('1.2720196495140689642524224617376+0.78615137775742328606955858584295892952312205783772323766490213j',(100,200))

Exceptions are normally raised in Python when the result of a real operation
is not defined over the reals; for example, ``sqrt(-4)`` will raise an
exception. The default context in gmpy2 implements the same behavior but by
setting allow_complex to True, complex results will be returned.

::

    >>> import gmpy2
    >>> from gmpy2 import mpc
    >>> gmpy2.sqrt(-4)
    mpfr('nan')
    >>> gmpy2.get_context(allow_complex=True)
    context(precision=53, real_prec=Default, imag_prec=Default,
            round=RoundToNearest, real_round=Default, imag_round=Default,
            emax=1073741823, emin=-1073741823,
            subnormalize=False,
            trap_underflow=False, underflow=False,
            trap_overflow=False, overflow=False,
            trap_inexact=False, inexact=False,
            trap_invalid=False, invalid=True,
            trap_erange=False, erange=False,
            trap_divzero=False, divzero=False,
            trap_expbound=False,
            allow_complex=True)
    >>> gmpy2.sqrt(-4)
    mpc('0.0+2.0j')

mpc Methods
-----------

**conjugate()**
    Returns the complex conjugate.

**digits()**
    Returns a two element tuple where each element represents the real and
    imaginary components as a 3-tuple containing the mantissa, the exponent,
    and the number of bits of precision. The mantissa is represented as a
    string in the specified base with up to 'prec' digits. If 'prec' is 0, as
    many digits that are available are returned. No more digits than available
    given x's precision are returned. 'base' must be between 2 and 62,
    inclusive.

mpc Attributes
--------------

**imag**
    Returns the imaginary component.

**precision**
    Returns a 2-tuple containing the precision of the real and imaginary
    components.

**rc**
    Returns a 2-tuple containing the ternary value of the real and imaginary
    components. The ternary value is 0 if the value of the component is exactly
    equal to the exact, infinite precision value. If the result code is 1, then
    the value of the component is greater than the exact value. If the result
    code is -1, then the value of the component is less than the exact,
    infinite precision value.

**real**
    Returns the real component.

mpc Functions
-------------

**acos(...)**
    acos(x) returns the arc-cosine of x.

**acosh(...)**
    acosh(x) returns the inverse hyperbolic cosine of x.

**add(...)**
    add(x, y) returns x + y. The type of the result is based on the types of
    the arguments.

**asin(...)**
    asin(x) returns the arc-sine of x.

**asinh(...)**
    asinh(x) return the inverse hyperbolic sine of x.

**atan(...)**
    atan(x) returns the arc-tangent of x.

**atanh(...)**
    atanh(x) returns the inverse hyperbolic tangent of x.

**cos(...)**
    cos(x) returns the cosine of x.

**cosh(...)**
    cosh(x) returns the hyperbolic cosine of x.

**div(...)**
    div(x, y) returns x / y. The type of the result is based on the types of
    the arguments.

**div_2exp(...)**
    div_2exp(x, n) returns an 'mpfr' or 'mpc' divided by 2**n.

**exp(...)**
    exp(x) returns e**x.

**fma(...)**
    fma(x, y, z) returns correctly rounded result of (x * y) + z.

**fms(...)**
    fms(x, y, z) returns correctly rounded result of (x * y) - z.

**is_inf(...)**
    is_inf(x) returns True if either the real or imaginary component of x is
    Infinity or -Infinity.

**is_nan(...)**
    is_nan(x) returns True if either the real or imaginary component of x is
    NaN (Not-A-Number).

**is_zero(...)**
    is_zero(x) returns True if x is zero.

**log(...)**
    log(x) returns the natural logarithm of x.

**log10(...)**
    log10(x) returns the base-10 logarithm of x.

**mpc(...)**
    mpc() returns an *mpc* object set to 0.0+0.0j.

    mpc(c[, precision=0]) returns a new 'mpc' object from an existing complex
    number (either a Python complex object or another 'mpc' object). If the
    precision is not specified, then the precision is taken from the current
    context. The rounding mode is always taken from the current context.

    mpc(r[, i=0[, precision=0]]) returns a new 'mpc' object by converting two
    non-complex numbers into the real and imaginary components of an 'mpc'
    object. If the precision is not specified, then the precision is taken from
    the current context. The rounding mode is always taken from the current
    context.

    mpc(s[, [precision=0[, base=10]]) returns a new 'mpc' object by converting
    a string s into a complex number. If base is omitted, then a base-10
    representation is assumed otherwise a base between 2 and 36 can be
    specified. If the precision is not specified, then the precision is taken
    from the current context. The rounding mode is always taken from the
    current context.

    In addition to the standard Python string representation of a complex
    number: ``"1+2j"``, the string representation used by the MPC library:
    ``"(1 2)"`` is also supported.

    .. note::
        The precision can be specified either a single number that is used for
        both the real and imaginary components, or as a 2-tuple that can
        specify different precisions for the real and imaginary components.

**mpc_random(...)**
    mpfc_random(random_state) returns a uniformly distributed number in the
    unit square [0,1]x[0,1]. The parameter *random_state* must be created by
    random_state() first.

**mul(...)**
    mul(x, y) returns x * y. The type of the result is based on the types of
    the arguments.

**mul_2exp(...)**
    mul_2exp(x, n) returns 'mpfr' or 'mpc' multiplied by 2**n.

**norm(...)**
    norm(x) returns the norm of a complex x. The norm(x) is defined as
    x.real**2 + x.imag**2. abs(x) is the square root of norm(x).

**phase(...)**
    phase(x) returns the phase angle, also known as argument, of a complex x.

**polar(...)**
    polar(x) returns the polar coordinate form of a complex x that is in
    rectangular form.

**proj(...)**
    proj(x) returns the projection of a complex x on to the Riemann sphere.

**rect(...)**
    rect(x) returns the polar coordinate form of a complex x that is in
    rectangular form.

**root_of_unity(...)**
    root_of_unity(n, k) returns the n-th root of mpc(1) raised to the k-th
    power. Requires MPC 1.1.0 or greater.

**sin(...)**
    sin(x) returns the sine of x.

**sinh(...)**
    sinh(x) returns the hyberbolic sine of x.

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

mpc Formatting
--------------

The *mpc* type supports the __format__() special method to allow custom output
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
    |     optional width.real_precision.imag_precision
    |     optional rounding mode:
    |        'U' -> round toward plus infinity
    |        'D' -> round toward minus infinity
    |        'Z' -> round toward zero
    |        'N' -> round to nearest
    |     optional output style:
    |        'P' -> Python style, 1+2j, (default)
    |        'M' -> MPC style, (1 2)
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







