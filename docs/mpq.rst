Multiple-precision Rationals
============================

gmpy2 provides a rational type call *mpq*. It should be a replacement for
Python's fractions.Fraction module.

::

    >>> import gmpy2
    >>> from gmpy2 import mpq
    >>> mpq(1,7)
    mpq(1,7)
    >>> mpq(1,7) * 11
    mpq(11,7)
    >>> mpq(11,7)/13
    mpq(11,91)

mpq Methods
-----------

**digits(...)**
    x.digits([base=10]) returns a Python string representing *x* in the
    given base (2 to 62, default is 10). A leading '-' is present if *x* < 0,
    but no leading '+' is present if *x* >= 0.

mpq Attributes
--------------

**denominator**
    x.denominator returns the denominator of *x*.

**numerator**
    x.numerator returns the numerator of *x*.

mpq Functions
-------------

**add(...)**
    add(x, y) returns *x* + *y*. The result type depends on the input
    types.

**div(...)**
    div(x, y) returns *x* / *y*. The result type depends on the input
    types.

**f2q(...)**
    f2q(x[, err]) returns the best *mpq* approximating *x* to within
    relative error *err*. Default is the precision of *x*. If *x* is not an
    *mpfr*, it is converted to an *mpfr*. Uses Stern-Brocot tree to find the
    best approximation. An *mpz* is returned if the denominator is 1. If
    *err* < 0, then the relative error sought is 2.0 ** *err*.

**mpq(...)**
    mpq() returns an *mpq* object set to 0/1.

    mpq(n) returns an *mpq* object with a numeric value *n*. Decimal and
    Fraction values are converted exactly.

    mpq(n, m) returns an *mpq* object with a numeric value *n* / *m*.

    mpq(s[, base=10]) returns an *mpq* object from a string *s* made up of
    digits in the given base. *s* may be made up of two numbers in the same
    base separated by a '/' character. If *base* == 10, then an embedded '.'
    indicates a number with a decimal fractional part.

**mul(...)**
    mul(x, y) returns *x* \* *y*. The result type depends on the input
    types.

**qdiv(...)**
    qdiv(x[, y=1]) returns *x/y* as *mpz* if possible, or as *mpq* if *x*
    is not exactly divisible by *y*.

**sub(...)**
    sub(x, y) returns *x* - *y*. The result type depends on the input
    types.




