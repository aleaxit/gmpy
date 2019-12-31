Overview of gmpy2
=================

Tutorial
--------

The *mpz* type is compatible with Python's built-in int/long type but is
significantly faster for large values. The cutover point for performance varies,
but can be as low as 20 to 40 digits. A variety of additional integer functions
are provided.
::

    >>> import gmpy2
    >>> from gmpy2 import mpz,mpq,mpfr,mpc
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

The *mpq* type is compatible with the *fractions.Fraction* type included with
Python.
::

    >>> mpq(3,7)/7
    mpq(3,49)
    >>> mpq(45,3) * mpq(11,8)
    mpq(165,8)

The most significant new features in gmpy2 are support for correctly rounded
arbitrary precision real and complex arithmetic based on the MPFR and MPC
libraries. Floating point contexts are used to control exceptional conditions.
For example, division by zero can either return an Infinity or raise an exception.
::

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
            trap_expbound=False,
            allow_complex=False)
    >>> mpfr(1)/0
    mpfr('inf')
    >>> gmpy2.get_context().trap_divzero=True
    >>> mpfr(1)/0
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    gmpy2.DivisionByZeroError: 'mpfr' division by zero in division
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
            trap_expbound=False,
            allow_complex=False)
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
    >>>


Miscellaneous gmpy2 Functions
-----------------------------

**from_binary(...)**
    from_binary(bytes) returns a gmpy2 object from a byte sequence created by
    to_binary().

**get_cache(...)**
    get_cache() returns the current cache size (number of objects) and the
    maximum size per object (number of limbs).

    gmpy2 maintains an internal list of freed *mpz*, *xmpz*, *mpq*, *mpfr*, and
    *mpc* objects for reuse. The cache significantly improves performance but
    also increases the memory footprint.

**license(...)**
    license() returns the gmpy2 license information.

**mp_limbsize(...)**
    mp_limbsize() returns the number of bits per limb used by the GMP or MPIR
    library.

**mp_version(...)**
    mp_version() returns the version of the GMP or MPIR library.

**mpc_version(...)**
    mpc_version() returns the version of the MPC library.

**mpfr_version(...)**
    mpfr_version() returns the version of the MPFR library.

**random_state(...)**
    random_state([seed]) returns a new object containing state information for
    the random number generator. An optional integer argument can be specified
    as the seed value. Only the Mersenne Twister random number generator is
    supported.

**set_cache(...)**
    set_cache(number, size) updates the maximum number of freed objects of each
    type that are cached and the maximum size (in limbs) of each object. The
    maximum number of objects of each type that can be cached is 1000. The
    maximum size of an object is 16384. The maximum size of an object is
    approximately 64K on 32-bit systems and 128K on 64-bit systems.

    .. note::
        The caching options are global to gmpy2. Changes are not thread-safe. A
        change in one thread will impact all threads.

**to_binary(...)**
    to_binary(x) returns a byte sequence from a gmpy2 object. All object types
    are supported.

**version(...)**
    version() returns the version of gmpy2.
