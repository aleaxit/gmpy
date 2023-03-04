Multiple-precision Rationals
============================

.. currentmodule:: gmpy2

gmpy2 provides a rational type `mpq`. It should be a replacement for
Python's `~fractions.Fraction` class.

.. doctest::

    >>> from gmpy2 import mpq
    >>> mpq(1,7)
    mpq(1,7)
    >>> mpq(1,7) * 11
    mpq(11,7)
    >>> mpq(11,7)/13
    mpq(11,91)

mpq type
--------

.. autoclass:: mpq

mpq Functions
-------------

.. autofunction:: qdiv
