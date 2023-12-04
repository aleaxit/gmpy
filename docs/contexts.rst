.. _contexts:

Contexts
========

.. currentmodule:: gmpy2

`context()` creates a new context.  `set_context()` will set the active
context.  `get_context()` will return a reference to the active context.  Note
that contexts are mutable: modifying the reference returned by `get_context()`
will modify the active context until a new context is enabled with
`set_context()`.  The `context.copy()` method will return a copy of the
context.  Contexts that implement the standard *single*, *double*, and
*quadruple* precision floating point types can be created using `ieee()`.

Context Type
------------

.. autoclass:: context

Context Functions
-----------------

.. autofunction:: get_context
.. autofunction:: ieee
.. autofunction:: local_context
.. autofunction:: set_context

Contexts and the with statement
-------------------------------

Contexts can also be used in conjunction with Python's :keyword:`with`
statement to temporarily change the context settings for a block of code and
then restore the original settings when the block of code exits.

`local_context` saves the current context and then creates a new context
based on the context passed as the first argument, or the current context
if no context is passed. The new context is modified if any optional keyword
arguments are given. The original active context is restored when the block
completes.

In the following example, the current context is saved by `local_context`
and then the block begins with a copy of the default context and the precision
set to 100. When the block is finished, the original context is restored.

.. doctest::

    >>> import gmpy2
    >>> print(gmpy2.sqrt(2))
    1.4142135623730951
    >>> with gmpy2.local_context(gmpy2.context(), precision=100):
    ...   print(gmpy2.sqrt(2))
    ...   gmpy2.get_context().precision += 100
    ...   print(gmpy2.sqrt(2))
    ...
    1.4142135623730950488016887242092
    1.4142135623730950488016887242096980785696718753769480731766796
    >>> print(gmpy2.sqrt(2))
    1.4142135623730951

Beginning in gmpy2 2.2.0, the :keyword:`with` statement can directly without
the use of `local_context`.

.. doctest::

    >>> import gmpy2
    >>> print(gmpy2.sqrt(2))
    1.4142135623730951
    >>> with gmpy2.context(precision=100):
    ...   print(gmpy2.sqrt(2))
    ...   gmpy2.get_context().precision += 100
    ...   print(gmpy2.sqrt(2))
    ...
    1.4142135623730950488016887242092
    1.4142135623730950488016887242096980785696718753769480731766796
    >>> print(gmpy2.sqrt(2))
    1.4142135623730951
