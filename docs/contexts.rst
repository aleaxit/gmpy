.. _contexts:

Contexts
========

.. currentmodule:: gmpy2

`context()` creates a new context with all options set to default.
`set_context()` will set the active context.  `get_context()` will
return a reference to the active context. Note that contexts are mutable:
modifying the reference returned by `get_context()` will modify the active
context until a new context is enabled with `set_context()`. The
`~context.copy()` method of a context will return a copy of the context.
Contexts that implement the standard *single* (32), *double* (64), and
*quadruple* (128) precision floating point types can be created using `ieee()`.

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

Nested contexts and coding style
--------------------------------

The ``with ... as ctx:`` coding style is commonly used. It does not behave
as expected with nested contexts. The first example shows the unexpected
behavior. The precision of sqrt(2) is correct since the global context
is properly undated. But ctx is not updated when the block exits. You
should always use get_context() to refer to the currently active context.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import ieee, sqrt, get_context
    >>> with ieee(64) as ctx:
    ...   print(ctx.precision, sqrt(2))
    ...   with ieee(128) as ctx:
    ...     print(ctx.precision, sqrt(2))
    ...   print(ctx.precision, sqrt(2))
    ...
    53 1.4142135623730951
    113 1.41421356237309504880168872420969798
    113 1.4142135623730951

This example uses the recommended coding style.

.. doctest::

    >>> import gmpy2
    >>> from gmpy2 import ieee, sqrt, get_context
    >>> with ieee(64):
    ...   print(get_context().precision, sqrt(2))
    ...   with ieee(128):
    ...     print(get_context().precision, sqrt(2))
    ...   print(get_context().precision, sqrt(2))
    ...
    53 1.4142135623730951
    113 1.41421356237309504880168872420969798
    53 1.4142135623730951
