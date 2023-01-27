Multiple-precision Rationals
============================

.. currentmodule:: gmpy2

gmpy2 provides a rational type :class:`mpq`. It should be a replacement for
Python's :class:`~fractions.Fraction` class.

::

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
   :members:

mpq Functions
-------------

.. autofunction:: add
   :noindex:
.. autofunction:: div
   :noindex:
.. autofunction:: f2q
.. autofunction:: mul
   :noindex:
.. autofunction:: qdiv
.. autofunction:: sub
   :noindex:
