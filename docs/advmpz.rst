Multiple-precision Integers (Advanced topics)
=============================================

The xmpz type
-------------

gmpy2 provides access to an experimental integer type called *xmpz*. The
*xmpz* type is a mutable integer type. In-place operations (+=, //=, etc.)
modify the original object and do not create a new object. Instances of
*xmpz* cannot be used as dictionary keys.

::

    >>> import gmpy2
    >>> from gmpy2 import xmpz
    >>> a = xmpz(123)
    >>> b = a
    >>> a += 1
    >>> a
    xmpz(124)
    >>> b
    xmpz(124)

The ability to change an *xmpz* object in-place allows for efficient and rapid
bit manipulation.

Individual bits can be set or cleared::

    >>> a[10]=1
    >>> a
    xmpz(1148)

Slice notation is supported. The bits referenced by a slice can be either 'read
from' or 'written to'. To clear a slice of bits, use a source value of 0. In
2s-complement format, 0 is represented by an arbitrary number of 0-bits. To set
a slice of bits, use a source value of ~0. The *tilde* operator inverts, or
complements the bits in an integer. (~0 is -1 so you can also use -1.) In
2s-complement format, -1 is represented by an arbitrary number of 1-bits.

If a value for *stop* is specified in a slice assignment and the actual
bit-length of the *xmpz* is less than *stop*, then the destination *xmpz* is
logically padded with 0-bits to length *stop*.

::

    >>> a=xmpz(0)
    >>> a[8:16] = ~0
    >>> bin(a)
    '0b1111111100000000'
    >>> a[4:12] = ~a[4:12]
    >>> bin(a)
    '0b1111000011110000'

Bits can be reversed::

    >>> bin(a)
    '0b10001111100'
    >>> a[::] = a[::-1]
    >>> bin(a)
    '0b111110001'

The *iter_bits()* method returns a generator that returns True or False for each
bit position. The methods *iter_clear()*, and *iter_set()* return generators
that return the bit positions that are 1 or 0. The methods support arguments
*start* and *stop* that define the beginning and ending bit positions that are
used. To mimic the behavior of slices. the bit positions checked include *start*
but the last position checked is *stop* - 1.

::

    >>> a=xmpz(117)
    >>> bin(a)
    '0b1110101'
    >>> list(a.iter_bits())
    [True, False, True, False, True, True, True]
    >>> list(a.iter_clear())
    [1, 3]
    >>> list(a.iter_set())
    [0, 2, 4, 5, 6]
    >>> list(a.iter_bits(stop=12))
    [True, False, True, False, True, True, True, False, False, False, False, False]

The following program uses the Sieve of Eratosthenes to generate a list of
prime numbers.

::

    from __future__ import print_function
    import time
    import gmpy2

    def sieve(limit=1000000):
        '''Returns a generator that yields the prime numbers up to limit.'''

	# Increment by 1 to account for the fact that slices  do not include
	# the last index value but we do want to include the last value for
	# calculating a list of primes.
	sieve_limit = gmpy2.isqrt(limit) + 1
	limit += 1

	# Mark bit positions 0 and 1 as not prime.
	bitmap = gmpy2.xmpz(3)

	# Process 2 separately. This allows us to use p+p for the step size
	# when sieving the remaining primes.
	bitmap[4 : limit : 2] = -1

	# Sieve the remaining primes.
	for p in bitmap.iter_clear(3, sieve_limit):
	    bitmap[p*p : limit : p+p] = -1

	return bitmap.iter_clear(2, limit)

    if __name__ == "__main__":
        start = time.time()
        result = list(sieve())
        print(time.time() - start)
        print(len(result))


Advanced Number Theory Functions
--------------------------------

The following functions are based on mpz_lucas.c and mpz_prp.c by David
Cleaver.

A good reference for probable prime testing is
http://www.pseudoprime.com/pseudo.html

**is_bpsw_prp(...)**
    is_bpsw_prp(n) will return True if *n* is a Baillie-Pomerance-Selfridge-Wagstaff
    probable prime. A BPSW probable prime passes the is_strong_prp() test with base
    2 and the is_selfridge_prp() test.

**is_euler_prp(...)**
    is_euler_prp(n,a) will return True if *n* is an Euler (also known as
    Solovay-Strassen) probable prime to the base *a*.

    | Assuming:
    |     gcd(n, a) == 1
    |     n is odd
    |
    | Then an Euler probable prime requires:
    |    a**((n-1)/2) == 1 (mod n)

**is_extra_strong_lucas_prp(...)**
    is_extra_strong_lucas_prp(n,p) will return True if *n* is an extra strong
    Lucas probable prime with parameters (p,1).

    | Assuming:
    |     n is odd
    |     D = p*p - 4, D != 0
    |     gcd(n, 2*D) == 1
    |     n = s*(2**r) + Jacobi(D,n), s odd
    |
    | Then an extra strong Lucas probable prime requires:
    |     lucasu(p,1,s) == 0 (mod n)
    |      or
    |     lucasv(p,1,s) == +/-2 (mod n)
    |      or
    |     lucasv(p,1,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_fermat_prp(...)**
    is_fermat_prp(n,a) will return True if *n* is a Fermat probable prime to the
    base a.

    | Assuming:
    |     gcd(n,a) == 1
    |
    | Then a Fermat probable prime requires:
    |     a**(n-1) == 1 (mod n)

**is_fibonacci_prp(...)**
    is_fibonacci_prp(n,p,q) will return True if *n* is a Fibonacci
    probable prime with parameters (p,q).

    | Assuming:
    |     n is odd
    |     p > 0, q = +/-1
    |     p*p - 4*q != 0
    |
    | Then a Fibonacci probable prime requires:
    |     lucasv(p,q,n) == p (mod n).

**is_lucas_prp(...)**
    is_lucas_prp(n,p,q) will return True if *n* is a Lucas probable prime with
    parameters (p,q).

    | Assuming:
    |     n is odd
    |     D = p*p - 4*q, D != 0
    |     gcd(n, 2*q*D) == 1
    |
    | Then a Lucas probable prime requires:
    |     lucasu(p,q,n - Jacobi(D,n)) == 0 (mod n)

**is_selfridge_prp(...)**
    is_selfridge_prp(n) will return True if *n* is a Lucas probable prime with
    Selfidge parameters (p,q). The Selfridge parameters are chosen by finding
    the first element D in the sequence {5, -7, 9, -11, 13, ...} such that
    Jacobi(D,n) == -1. Let p=1 and q = (1-D)/4 and then perform a Lucas
    probable prime test.

**is_strong_bpsw_prp(...)**
    is_strong_bpsw_prp(n) will return True if *n* is a strong
    Baillie-Pomerance-Selfridge-Wagstaff probable prime. A strong BPSW
    probable prime passes the is_strong_prp() test with base 2 and the
    is_strongselfridge_prp() test.

**is_strong_lucas_prp(...)**
    is_strong_lucas_prp(n,p,q) will return True if *n* is a strong Lucas
    probable prime with parameters (p,q).

    | Assuming:
    |     n is odd
    |     D = p*p - 4*q, D != 0
    |     gcd(n, 2*q*D) == 1
    |     n = s*(2**r) + Jacobi(D,n), s odd
    |
    | Then a strong Lucas probable prime requires:
    |     lucasu(p,q,s) == 0 (mod n)
    |      or
    |     lucasv(p,q,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_strong_prp(...)**
    is_strong_prp(n,a) will return True if *n* is a strong (also known as
    Miller-Rabin) probable prime to the base a.

    | Assuming:
    |     gcd(n,a) == 1
    |     n is odd
    |     n = s*(2**r) + 1, with s odd
    |
    | Then a strong probable prime requires one of the following is true:
    |     a**s == 1 (mod n)
    |      or
    |     a**(s*(2**t)) == -1 (mod n) for some t, 0 <= t < r.

**is_strong_selfridge_prp(...)**
    is_strong_selfridge_prp(n) will return True if *n* is a strong Lucas
    probable prime with Selfidge parameters (p,q). The Selfridge parameters are
    chosen by finding the first element D in the sequence
    {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) == -1. Let p=1 and
    q = (1-D)/4 and then perform a strong Lucas probable prime test.

**lucasu(...)**
    lucasu(p,q,k) will return the k-th element of the Lucas U sequence defined
    by p,q. p*p - 4*q must not equal 0; k must be greater than or equal to 0.

**lucasu_mod(...)**
    lucasu_mod(p,q,k,n) will return the k-th element of the Lucas U sequence
    defined by p,q (mod n). p*p - 4*q must not equal 0; k must be greater than
    or equal to 0; n must be greater than 0.

**lucasv(...)**
    lucasv(p,q,k) will return the k-th element of the Lucas V sequence defined
    by parameters (p,q). p*p - 4*q must not equal 0; k must be greater than or
    equal to 0.

**lucasv_mod(...)**
    lucasv_mod(p,q,k,n) will return the k-th element of the Lucas V sequence
    defined by parameters (p,q) (mod n). p*p - 4*q must not equal 0; k must be
    greater than or equal to 0; n must be greater than 0.

