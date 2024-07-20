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
*quadruple* precision floating-point types can be created using `ieee()`.

Context Type
------------

.. autoclass:: context

Context Functions
-----------------

.. autofunction:: get_context
.. autofunction:: ieee
.. autofunction:: local_context
.. autofunction:: set_context
