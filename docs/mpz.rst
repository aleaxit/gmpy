Multiple-precision Integers
===========================

The gmpy2 *mpz* type supports arbitrary precision integers. It should be a
drop-in replacement for Python's *long* type. Depending on the platform and the
specific operation, an *mpz* will be faster than Python's *long* once the
precision exceeds 20 to 50 digits. All the special integer functions in GMP are
supported.

Examples
--------

::

    >>> import gmpy2
    >>> from gmpy2 import mpz
    >>> mpz('123') + 1
    mpz(124)
    >>> 10 - mpz(1)
    mpz(9)
    >>> gmpy2.is_prime(17)
    True
    >>> mpz('1_2')
    mpz(12)

.. note::
    The use of ``from gmpy2 import *`` is not recommended. The names in gmpy2
    have been chosen to avoid conflict with Python's builtin names but gmpy2
    does use names that may conflict with other modules or variable names.

.. note::
   gmpy2.mpz() ignores all embedded underscore characters. It does not attempt
   to be 100% compatible with all Python exceptions.

mpz Methods
-----------

**bit_clear(...)**
    x.bit_clear(n) returns a copy of *x* with bit *n* set to 0.

**bit_flip(...)**
    x.bit_flip(n) returns a copy of *x* with bit *n* inverted.

**bit_length(...)**
    x.bit_length() returns the number of significant bits in the radix-2
    representation of *x*. For compatibility with Python, mpz(0).bit_length()
    returns 0.

**bit_scan0(...)**
    x.bit_scan0(n) returns the index of the first 0-bit of *x* with
    index >= *n*. If there are no more 0-bits in *x* at or above index *n*
    (which can only happen for *x* < 0, assuming an infinitely long 2's
    complement format), then None is returned. *n* must be >= 0.

**bit_scan1(...)**
    x.bit_scan1(n) returns the index of the first 1-bit of *x* with
    index >= *n*. If there are no more 1-bits in *x* at or above index *n*
    (which can only happen for *x* >= 0, assuming an infinitely long 2's
    complement format), then None is returned. *n* must be >= 0.

**bit_set(...)**
    x.bit_set(n) returns a copy of *x* with bit *n* set to 1.

**bit_test(...)**
    x.bit_test(n) returns True if bit *n* of *x* is set, and False if it
    is not set.

**conjugtae(...)**
    Return the conjugate of x (which is just a new reference to x since x
    not a complex number).

**denominator(...)**
    x.denominator() returns mpz(1).

**digits(...)**
    x.digits([base=10]) returns a string representing *x* in radix *base*.

**imag**
    Return the imaginary component of an *mpz*. Always *mpz(0)*.

**is_congruent(...)**
    x.is_congruent(y, m) returns True if *x* is congruent to *y* modulo *m*,
    else returns False.

**is_divisible(...)**
    x.is_divisible(d) returns True if *x* is divisible by *d*, else returns
    False.

**is_even(...)**
    x.is_even() returns True if *x* is even, else returns False.

**is_odd(...)**
    x.is_odd() returns True if *x* is even, else returns False.

**is_power(...)**
    x.is_power() returns True if *x* is a perfect power (there exists integers
    *y* and *n* > 1, such that x=y**n), else returns False.

**is_prime(...)**
    x.is_prime() returns True if *x* is _probably_ prime, else False if *x* is
    definitely composite.

    See the documentation for *gmpy2.is_prime* for details on the underlaying
    primality tests that are performed.

**is_square(...)**
    x.is_square() returns True if *x* is a perfect square, else returns False.

**num_digits(...)**
    x.num_digits([base=10]) returns the length of the string representing
    the absolute value of *x* in radix *base*. The result is correct if base is
    a power of 2. For other bases, the result is usually correct but may
    be 1 too large. *base* can range between 2 and 62, inclusive.

**numerator(...)**
    x.numerator() returns a copy of *x*.

**real(...)**
    x.real returns a copy of *x*.

mpz Functions
-------------

**add(...)**
    add(x, y) returns *x* + *y*. The result type depends on the input
    types.

**bincoef(...)**
    bincoef(x, n) returns the binomial coefficient. *n* must be >= 0.

**bit_clear(...)**
    bit_clear(x, n) returns a copy of *x* with bit *n* set to 0.

**bit_flip(...)**
    bit_flip(x, n) returns a copy of *x* with bit *n* inverted.

**bit_length(...)**
    bit_length(x) returns the number of significant bits in the radix-2
    representation of *x*. For compatibility with Python, mpz(0).bit_length()
    returns 0 while mpz(0).num_digits(2) returns 1.

**bit_mask(...)**
    bit_mask(n) returns an *mpz* object exactly *n* bits in length with all
    bits set.

**bit_scan0(...)**
    bit_scan0(x, n) returns the index of the first 0-bit of *x* with
    index >= *n*. If there are no more 0-bits in *x* at or above index *n*
    (which can only happen for *x* < 0, assuming an infinitely long 2's
    complement format), then None is returned. *n* must be >= 0.

**bit_scan1(...)**
    bit_scan1(x, n) returns the index of the first 1-bit of *x* with
    index >= *n*. If there are no more 1-bits in *x* at or above index *n*
    (which can only happen for *x* >= 0, assuming an infinitely long 2's
    complement format), then None is returned. *n* must be >= 0.

**bit_set(...)**
    bit_set(x, n) returns a copy of *x* with bit *n* set to 1.

**bit_test(...)**
    bit_test(x, n) returns True if bit *n* of *x* is set, and False if it
    is not set.

**c_div(...)**
    c_div(x, y) returns the quotient of *x* divided by *y*. The quotient is
    rounded towards +Inf (ceiling rounding). *x* and *y* must be integers.

**c_div_2exp(...)**
    c_div_2exp(x, n) returns the quotient of *x* divided by 2**n. The
    quotient is rounded towards +Inf (ceiling rounding). *x* must be an integer
    and *n* must be > 0.

**c_divmod(...)**
    c_divmod(x, y) returns the quotient and remainder of *x* divided by
    *y*. The quotient is rounded towards +Inf (ceiling rounding) and the
    remainder will have the opposite sign of *y*. *x* and *y* must be integers.

**c_divmod_2exp(...)**
    c_divmod_2exp(x ,n) returns the quotient and remainder of *x* divided
    by 2**n. The quotient is rounded towards +Inf (ceiling rounding) and the
    remainder will be negative or zero. *x* must be an integer and *n* must
    be > 0.

**c_mod(...)**
    c_mod(x, y) returns the remainder of *x* divided by *y*. The remainder
    will have the opposite sign of *y*. *x* and *y* must be integers.

**c_mod_2exp(...)**
    c_mod_2exp(x, n) returns the remainder of *x* divided by 2**n. The
    remainder will be negative. *x* must be an integer and *n* must be > 0.

**comb(...)**
    comb(x, n) returns the number of combinations of *x* things, taking *n*
    at a time. *n* must be >= 0.

**digits(...)**
    digits(x[, base=10]) returns a string representing *x* in radix *base*.

**div(...)**
    div(x, y) returns *x* / *y*. The result type depends on the input
    types.

**divexact(...)**
    divexact(x, y) returns the quotient of *x* divided by *y*. Faster than
    standard division but requires the remainder is zero!

**divm(...)**
    divm(a, b, m) returns *x* such that *b* * *x* == *a* modulo *m*. Raises
    a ZeroDivisionError exception if no such value *x* exists.

**double_fac(...)**
    double_fac(n) returns the exact double factorial of *n*.

**f_div(...)**
    f_div(x, y) returns the quotient of *x* divided by *y*. The quotient
    is rounded towards -Inf (floor rounding). *x* and *y* must be integers.

**f_div_2exp(...)**
    f_div_2exp(x, n) returns the quotient of *x* divided by 2**n. The
    quotient is rounded towards -Inf (floor rounding). *x* must be an integer
    and *n* must be > 0.

**f_divmod(...)**
    f_divmod(x, y) returns the quotient and remainder of *x* divided by
    *y*. The quotient is rounded towards -Inf (floor rounding) and the
    remainder will have the same sign as *y*. *x* and *y* must be integers.

**f_divmod_2exp(...)**
    f_divmod_2exp(x, n) returns quotient and remainder after dividing *x*
    by 2**n. The quotient is rounded towards -Inf (floor rounding) and the
    remainder will be positive. *x* must be an integer and *n* must be > 0.

**f_mod(...)**
    f_mod(x, y) returns the remainder of *x* divided by *y*. The remainder
    will have the same sign as *y*. *x* and *y* must be integers.

**f_mod_2exp(...)**
    f_mod_2exp(x, n) returns remainder of *x* divided by 2**n. The
    remainder will be positive. *x* must be an integer and *n* must be > 0.

**fac(...)**
    fac(n) returns the exact factorial of *n*. Use factorial() to get the
    floating-point approximation.

**fib(...)**
    fib(n) returns the *n*-th Fibonacci number.

**fib2(...)**
    fib2(n) returns a 2-tuple with the (*n*-1)-th and *n*-th Fibonacci
    numbers.

**gcd(...)**
    gcd(a, b) returns the greatest common divisor of integers *a* and
    *b*.

**gcdext(...)**
    gcdext(a, b) returns a 3-element tuple (*g*, *s*, *t*) such that

    *g* == gcd(*a*, *b*) and *g* == *a* * *s*  + *b* * *t*

**hamdist(...)**
    hamdist(x, y) returns the Hamming distance (number of bit-positions
    where the bits differ) between integers *x* and *y*.

**invert(...)**
    invert(x, m) returns *y* such that *x* * *y* == 1 modulo *m*, or 0
    if no such *y* exists.

**iroot(...)**
    iroot(x,n) returns a 2-element tuple (*y*, *b*) such that *y* is the integer
    *n*-th root of *x* and *b* is True if the root is exact. *x* must be >= 0
    and *n* must be > 0.

**iroot_rem(...)**
    iroot_rem(x,n) returns a 2-element tuple (*y*, *r*) such that *y* is
    the integer *n*-th root of *x* and *x* = y**n + *r*. *x* must be >= 0 and
    *n* must be > 0.

**is_even(...)**
    is_even(x) returns True if *x* is even, False otherwise.

**is_odd(...)**
    is_odd(x) returns True if *x* is odd, False otherwise.

**is_power(...)**
    is_power(x) returns True if *x* is a perfect power, False otherwise.

**is_prime(...)**
    is_prime(x[, n=25]) returns True if *x* is **probably** prime. False
    is returned if *x* is definitely composite. *x* is checked for small
    divisors and up to *n* Miller-Rabin tests are performed. The actual tests
    performed may vary based on version of GMP or MPIR used.

**is_square(...)**
    is_square(x) returns True if *x* is a perfect square, False otherwise.

**isqrt(...)**
    isqrt(x) returns the integer square root of an integer *x*. *x* must be
    >= 0.

**isqrt_rem(...)**
    isqrt_rem(x) returns a 2-tuple (*s*, *t*) such that *s* = isqrt(*x*)
    and *t* = *x* - *s* * *s*. *x* must be >= 0.

**jacobi(...)**
    jacobi(x, y) returns the Jacobi symbol (*x* | *y*). *y* must be odd and
    > 0.

**kronecker(...)**
    kronecker(x, y) returns the Kronecker-Jacobi symbol (*x* | *y*).

**lcm(...)**
    lcm(a, b) returns the lowest common multiple of integers *a* and *b*.

**legendre(...)**
    legendre(x, y) returns the Legendre symbol (*x* | *y*). *y* is assumed
    to be an odd prime.

**lucas(...)**
    lucas(n) returns the *n*-th Lucas number.

**lucas2(...)**
    lucas2(n) returns a 2-tuple with the (*n*-1)-th and *n*-th Lucas
    numbers.

**mpz(...)**
    mpz() returns a new *mpz* object set to 0.

    mpz(n) returns a new *mpz* object from a numeric value *n*. If *n* is
    not an integer, it will be truncated to an integer.

    mpz(s[, base=0]) returns a new *mpz* object from a string *s* made of
    digits in the given base. If base = 0, then binary, octal, or hex Python
    strings are recognized by leading 0b, 0o, or 0x characters. Otherwise the
    string is assumed to be decimal. Values for base can range between 2 and 62.

**mpz_random(...)**
    mpz_random(random_state, n) returns a uniformly distributed random
    integer between 0 and *n*-1. The parameter *random_state* must be created
    by random_state() first.

**mpz_rrandomb(...)**
    mpz_rrandomb(random_state, b) returns a random integer between 0 and
    2**b - 1 with long sequences of zeros and one in its binary representation.
    The parameter *random_state* must be created by random_state() first.

**mpz_urandomb(...)**
    mpz_urandomb(random_state, b) returns a uniformly distributed random
    integer between 0 and 2**b - 1. The parameter *random_state* must be
    created by random_state() first.

**mul(...)**
    mul(x, y) returns *x* \* *y*. The result type depends on the input
    types.

**multi_fac(...)**
    multi_fac(n, m) returns the m-multi-factorial of *n* i.e n!^m.

**next_prime(...)**
    next_prime(x) returns the next **probable** prime number > *x*.

**num_digits(...)**
    num_digits(x[, base=10]) returns the length of the string representing
    the absolute value of *x* in radix *base*. The result is correct if base is
    a power of 2. For other bases, the result is usually correct but may
    be 1 too large. *base* can range between 2 and 62, inclusive.

**popcount(...)**
    popcount(x) returns the number of bits with value 1 in *x*. If *x* < 0,
    the number of bits with value 1 is infinite so -1 is returned in that case.

**powmod(...)**
    powmod(x, y, m) returns (*x* ** *y*) mod *m*. The exponent *y* can be
    negative, and the correct result will be returned if the inverse of *x*
    mod *m* exists. Otherwise, a ValueError is raised.

**primorial(...)**
    primorial(n) returns the exact primorial of *n*, i.e. the product of all
    positive prime numbers <= *n*.

**remove(...)**
    remove(x, f) will remove the factor *f* from *x* as many times as possible
    and return a 2-tuple (*y*, *m*) where *y* = *x* // (*f* ** *m*). *f* does
    not divide *y*. *m* is the multiplicity of the factor *f* in *x*. *f* must
    be > 1.

**sub(...)**
    sub(x, y) returns *x* - *y*. The result type depends on the input
    types.

**t_div(...)**
    t_div(x, y) returns the quotient of *x* divided by *y*. The quotient
    is rounded towards zero (truncation). *x* and *y* must be integers.

**t_div_2exp(...)**
    t_div_2exp(x, n) returns the quotient of *x* divided by 2**n. The
    quotient is rounded towards zero (truncation). *n* must be > 0.

**t_divmod(...)**
    t_divmod(x, y) returns the quotient and remainder of *x* divided by
    *y*. The quotient is rounded towards zero (truncation) and the remainder
    will have the same sign as *x*. *x* and *y* must be integers.

**t_divmod_2exp(...)**
    t_divmod_2exp(x, n) returns the quotient and remainder of *x* divided
    by 2**n. The quotient is rounded towards zero (truncation) and the
    remainder will have the same sign as *x*. *x* must be an integer and *n*
    must be > 0.

**t_mod(...)**
    t_mod(x, y) returns the remainder of *x* divided by *y*. The remainder
    will have the same sign as *x*. *x* and *y* must be integers.

**t_mod_2exp(...)**
    t_mod_2exp(x, n) returns the remainder of *x* divided by 2**n. The
    remainder will have the same sign as *x*. *x* must be an integer and *n*
    must be > 0.


