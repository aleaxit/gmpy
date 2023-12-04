from .gmpy2 import *
from .gmpy2 import __version__
# Internal variables/functions are not imported by * above.
# These are used by some python level functions and are needed
# at the top level.
# Use try...except to for static builds were _C_API is not available.
try:
    from .gmpy2 import _C_API, _mpmath_normalize, _mpmath_create
except ImportError:
    from .gmpy2 import _mpmath_normalize, _mpmath_create


def local_context(*args, **kwargs):
    """
    local_context(**kwargs) -> context
    local_context(context, /, **kwargs) -> context

    Create a context manager object that will restore the current context
    when the 'with ...' block terminates. The temporary context for the
    'with ...' block is based on the current context if no context is
    specified. Keyword arguments are supported and will modify the
    temporary new context.
    """
    import warnings
    warnings.warn("local_context() is deprecated, use context() instead",
                  DeprecationWarning)
    return context(*args, **kwargs)
