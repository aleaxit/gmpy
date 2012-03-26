Multiple-precision Integers (Advanced topics)
=============================================

The xmpz type
-------------

gmpy2 provides access to an experimental integer type called *xmpz*. The
*xmpz* type is a mutable integer type. Instances of *xmpz* cannot be used as
dictionary keys. In-place operations (+=, //=, etc.) modify the orignal object
and do not create a new object. The ability to change an *xmpz* object in-place
allows for efficient and rapid bit manipulation.

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

Individual bits can be set or cleared::

    >>> a[10]=1
    >>> a
    xmpz(1148)

Bits can be reversed::

    >>> bin(a)
    '0b10001111100'
    >>> a[::] = a[::-1]
    >>> bin(a)
    '0b111110001'

Advanced Number Theory Functions
--------------------------------

The following functions are based on mpz_lucas.c and mpz_prp.c by David
Cleaver.

A good reference for pseudo-prime testing is
http://www.pseudoprime.com/pseudo.html

**is_bpsw_prp(...)**
    is_bpsw_prp(n) will return True if *n* is a Baillie-Pomerance-Selfridge-Wagstaff
    pseudo-prime. A BPSW pseudoprime passes the is_strong_prp() test with base
    2 and the is_selfridge_prp() test.

**is_euler_prp(...)**
    is_euler_prp(n,a) will return True if *n* is an Euler (also known as
    Solovay-Strassen) pseudo-prime to the base *a*.

    | Assuming:
    |     gcd(n, a) == 1
    |     n is odd
    |
    | Then an Euler pseudo-prime requires:
    |    a**((n-1)/2) == 1 (mod n)

**is_extra_strong_lucas_prp(...)**
    is_extra_strong_lucas_prp(n,p) will return True if *n* is an extra strong
    Lucas pseudo-prime with parameters (p,q).

    | Assuming:
    |     n is odd
    |     D = p*p - 4, D != 0
    |     gcd(n, 2*D) == 1
    |     n = s*(2**r) + Jacobi(D,n), s odd
    |
    | Then an extra strong Lucas pseudoprime requires:
    |     lucasu(p,1,s) == 0 (mod n)
    |      or
    |     lucasv(p,1,s) == +/-2 (mod n)
    |      or
    |     lucasv(p,1,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_fermat_prp(...)**
    is_fermat_prp(n,a) will return True if *n* is a Fermat pseudo-prime to the
    base a.

    | Assuming:
    |     gcd(n,a) == 1
    |
    | Then a Fermat pseudoprime requires:
    |     a**(n-1) == 1 (mod n)

**is_fibonacci_prp(...)**
    is_fibonacci_prp(n,p,q) will return True if *n* is an Fibonacci
    pseudo-prime with parameters (p,q).

    | Assuming:
    |     n is odd
    |     p > 0, q = +/-1
    |     p*p - 4*q != 0
    |
    | Then a Fibonacci pseudo-prime requires:
    |     lucasv(p,q,n) == p (mod n).

**is_lucas_prp(...)**
    is_lucas_prp(n,p,q) will return True if *n* is a Lucas pseudo-prime with
    parameters (p,q).

    | Assuming:
    |     n is odd
    |     D = p*p - 4*q, D != 0
    |     gcd(n, 2*q*D) == 1
    |
    | Then a Lucas pseudo-prime requires:
    |     lucasu(p,q,n - Jacobi(D,n)) == 0 (mod n)

**is_selfridge_prp(...)**
    is_selfridge_prp(n) will return True if *n* is a Lucas pseudo-prime with
    Selfidge parameters (p,q). The Selfridge parameters are chosen by finding
    the first element D in the sequence {5, -7, 9, -11, 13, ...} such that
    Jacobi(D,n) == -1. Let p=1 and q = (1-D)/4 and then perform a Lucas
    pseudo-prime test.

**is_strong_bpsw_prp(...)**
    is_strong_bpsw_prp(n) will return True if *n* is a strong
    Baillie-Pomerance-Selfridge-Wagstaff pseudo-prime. A strong BPSW
    pseudo-prime passes the is_strong_prp() test with base 2 and the
    is_strongselfridge_prp() test.

**is_strong_lucas_prp(...)**
    is_strong_lucas_prp(n,p,q) will return True if *n* is a strong Lucas
    pseudo-prime with parameters (p,q).

    | Assuming:
    |     n is odd
    |     D = p*p - 4*q, D != 0
    |     gcd(n, 2*q*D) == 1
    |     n = s*(2**r) + Jacobi(D,n), s odd
    |
    | Then a strong Lucas pseudoprime requires:
    |     lucasu(p,q,s) == 0 (mod n)
    |      or
    |     lucasv(p,q,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r

**is_strong_prp(...)**
    is_strong_prp(n,a) will return True if *n* is an strong (also known as
    Miller-Rabin) pseudo-prime to the base a.

    | Assuming:
    |     gcd(n,a) == 1
    |     n is odd
    |     n = s*(2**r) + 1, with s odd
    |
    | Then a strong pseudoprime requires one of the following is true:
    |     a**s == 1 (mod n)
    |      or
    |     a**(s*(2**t)) == -1 (mod n) for some t, 0 <= t < r.

**is_strong_selfridge_prp(...)**
    is_strong_selfridge_prp(n) will return True if *n* is a strong Lucas
    pseudo-prime with Selfidge parameters (p,q). The Selfridge parameters are
    chosen by finding the first element D in the sequence
    {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) == -1. Let p=1 and
    q = (1-D)/4 and then perform a strong Lucas pseudo-prime test.

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

