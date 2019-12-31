Conversion methods and gmpy2's numbers
======================================

Conversion methods
------------------

A python object could interact with gmpy2 if it implements one of the following methods:

- **__mpz__** : return an object of <type 'mpz'>.
- **__mpq__** : return an object of <type 'mpq'>.
- **__mpfr__** : return an object of <type 'mpfr'>.
- **__mpc__** : return an object of <type 'mpc'>.

| Implementing on of these methods allow gmpy2 to convert a python object into a gmpy2 type.
| Example::

    >>> from gmpy2 import mpz
    >>> class CustInt:
    ...     def __init__(self, x):
    ...             self.x = x
    ...     def __mpz__(self):
    ...             return mpz(self.x)
    ...
    >>> ci = CustInt(5)
    >>> z = mpz(ci); z
    mpz(5)
    >>> type(z)
    <type 'mpz'>

Arithmetic operations
---------------------

| gmpy2 allow arithmetic operations between gmpy2 numbers and objects with conversion methods.
| Operation with object that implements floating conversion and exact conversion methods are not supported.
| That means that only the following cases are supported:

- An integer type have to implement **__mpz__**
- A rational type have to implement **__mpq__** and can implement **__mpz__**
- A real type have to implement **__mpfr__**
- A complex type have to implement **__mpc__** and can implement **__mpfr__**

Examples::

    >>> from gmpy2 import mpz, mpq, mpfr, mpc
    >>> class Q:
    ...     def __mpz__(self): return mpz(1)
    ...     def __mpq__(self): return mpq(3,2)
    >>> q = Q()
    >>> mpz(2) + q
    mpq(7,2)
    >>> mpq(1,2) * q
    mpq(3,4)
    >>> mpfr(10) * q
    mpfr('15.0')
