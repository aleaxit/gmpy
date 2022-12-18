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

**bit_count(...)**
    x.bit_count() returns the number of 1-bits set in abs(x).

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

**bit_count(...)**
    bit_count(x) returns a the number of 1 bits in the binary 
    representation of *x*. Differs from popcount() for x <0.

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
    gcd(...) returns the greatest common multiple of a sequence of integers.

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

**is_bpsw_prp(...)**
    is_bpsw_prp(n) returns True if *n* is a Baillie-Pomerance-Selfridge-Wagstaff
    probable prime. A BPSW probable prime passes both the is_strong_prp() test 
    with base 2 and the is_selfridge_prp() test.

**is_congruent(...)**
    is_congurent(x, y, m) returns True if *x* is congruent to *y* modulo *m*, 
    else return False.

**is_divisible(...)**
    is_divisible(x, d) returns True if *x* is divisible by *d*, else return False.

**is_euler_prp(...)**
    is_euler_prp(n, a) returns True if *n* is an Euler probable prime to the
    base *a*.

    Assuming:
        gcd(n, a) == 1
        n is odd
    then "n* is an Euler prp if:
        a**((n-1)/2) == 1 (mod n)

**is_even(...)**
    is_even(x) returns True if *x* is even, False otherwise.

**is_extra_strong_lucas_prp(...)**
    is_extra_strong_lucas_prp(n, p) returns True if n is an extra strong Lucas probable
    prime with parameters (p, 1). 
    
    Assuming:
        n is odd
        D = p*p - 4
        D != 0
        gcd(n, 2*D) == 1
        n = s*(2**r) + Jacobi(D,n), s odd
    then *n* is an extra strong Lucas probable prime if:
        lucasu(p,1,s) == 0 (mod n)
        or
        lucasv(p,1,s) == +/-2 (mod n)
        or
        lucasv(p,1,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_fermat_prp(...)**
    is_fermat_prp(n ,a) returns True if *n* is a Fermat probable prime
    to the base *a*.

    Assuming:
        gcd(n,a) == 1
    then *n* is a Fermat probable prime if:
        a**(n-1) == 1 (mod n)

**is_fibonacci_prp(...)**
    is_fibonacci_prp(n,p,q) returns True if *n* is a Fibonacci probable prime
    with parameters (p,q).

    Assuming:
        n is odd
        p > 0
        q = +/-1
        p*p - 4*q != 0
    then *n* is a Fibonacci probable prime if:
        lucasv(p,q,n) == p (mod n).

**is_lucas_prp(...)**
    is_lucas_prp(n,p,q) returns True if *n* is a Lucas probable prime with 
    parameters (p,q).

    Assuming:
        n is odd
        D = p*p - 4*q, D != 0
        gcd(n, 2*q*D) == 1
    then *n* is a Lucas probable prime if:
        lucasu(p,q,n - Jacobi(D,n)) == 0 (mod n)

**is_odd(...)**
    is_odd(x) returns True if *x* is odd, False otherwise.

**is_power(...)**
    is_power(x) returns True if *x* is a perfect power, False otherwise.

**is_prime(...)**
    is_prime(x[, n=25]) returns True if *x* is **probably** prime. False
    is returned if *x* is definitely composite. *x* is checked for small
    divisors and up to *n* Miller-Rabin tests are performed. The actual tests
    performed may vary based on version of GMP used.

**is_selfridge_prp(...)**
    is_selfridge_prp(n) returns True if *n* is a Lucas probable prime with
    Selfidge parameters (p,q). The Selfridge parameters are chosen by finding
    the first element D in the sequence {5, -7, 9, -11, 13, ...} such that
    Jacobi(D,n) == -1. Then let p=1 and q=(1-D)/4 and perform the Lucas probable
    prime test.

**is_square(...)**
    is_square(x) returns True if *x* is a perfect square, False otherwise.

**is_strong_bpsw_prp(...)**
    is_strong_bpsw_prp(n) returns True if *n* is a strong Baillie-Pomerance-
    Selfridge-Wagstaff probable prime. A strong BPSW probable prime passes the
    is_strong_prp() test with base 2 and the is_strong_selfridge_prp() test.

**is_strong_lucas_prp(...)**
    is_strong_lucas_prp(n,p,q) returns True if *n* is a strong Lucas probable
    prime with parameters (p,q).

    Assuming:
        n is odd
        D = p*p - 4*q, D != 0
        gcd(n, 2*q*D) == 1
        n = s*(2**r) + Jacobi(D,n), s odd
    then *n* is a strong Lucas probable prime if:
        lucasu(p,q,s) == 0 (mod n)
        or
        lucasv(p,q,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_strong_prp(...)**
    is_strong_prp(n,a) returns True if *n* is a strong (also known as
    Miller-Rabin) probable prime to the base *a*.

    Assuming:
        gcd(n,a) == 1
        n is odd
        n = s*(2**r) + 1, with s odd
    then *n* is a strong probable prime if:
        a**s == 1 (mod n)
        or
        a**(s*(2**t)) == -1 (mod n) for some t, 0 <= t < r.

**is_strong_selfridge_prp(...)**
    is_strong_selfridge_prp(n) returns True if *n* is a strong Lucas
    probable prime with Selfidge parameters (p,q). The Selfridge parameters
    are chosen by finding the first element D in the sequence {5, -7, 9, -11,
    13, ...} such that Jacobi(D,n) == -1. Then let p=1 and q = (1-D)/4 and
    perform a strong Lucas probable prime test.

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
    lcm(...) returns the lowest common multiple of a sequence of integers.

**legendre(...)**
    legendre(x, y) returns the Legendre symbol (*x* | *y*). *y* is assumed
    to be an odd prime.

**lucas(...)**
    lucas(n) returns the *n*-th Lucas number.

**lucas2(...)**
    lucas2(n) returns a 2-tuple with the (*n*-1)-th and *n*-th Lucas
    numbers.

**lucasu(...)**
    lucasu(p,q,k) returns the *k*-th element of the Lucas U sequence defined
    by (p,q). p*p - 4*q must not equal 0; *k* must be greater than or equal to 0.

**lucasu_mod(...)**
    lucasu_mod(p,q,k,n) returns the *k*-th element of the Lucas U sequence defined
    by (p,q) modulo *n*. p*p - 4*q must not equal 0; *k* must be greater than or
    equal to 0; *n* must be greater than 0.

**lucasv(...)**
    lucasv(p,q,k) returns the *k*-th element of the Lucas V sequence defined
    by (p,q). p*p - 4*q must not equal 0; *k* must be greater than or equal to 0.

**lucasv_mod(...)**
    lucasv_mod(p,q,k,n) returns the *k*-th element of the Lucas V sequence defined
    by (p,q) modulo *n*. p*p - 4*q must not equal 0; *k* must be greater than or
    equal to 0; *n* must be greater than 0.

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

**powmod_exp_list(...)**
    powmod_exp_list(base, exp_lst, mod) returns list(powmod(base, i, mod) for i in exp_lst).
    Releases the GIL so can be easily run in multiple threads. 
    
    Experimental in gmpy2 2.1.x. The capability will continue to exist in future
    versions but the name may change.

**powmod_base_list(...)**
    powmod_base_list(base_list, exp, mod) returns list(powmod(i, exp, mod) for i in lst).
    Releases the GIL so can be easily run in multiple threads. 
    
    Experimental in gmpy2 2.1.x. The capability will continue to exist in future
    versions but the name may change.

**powmod_sec(...)**
    powmod_sec(x, y, m) returns (*x* ** *y*) mod *m*. The calculation uses a
    constant time algorithm to reduce the risk of side channel attacks. *y* must
    be an integer >0. *m* must be an odd integer.

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
