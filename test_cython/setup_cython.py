from setuptools import Extension, setup
from Cython.Build import cythonize
import platform
import sys
import os
import gmpy2

ON_WINDOWS = platform.system() == 'Windows'

extensions = [
    Extension("test_cython", ["test_cython.pyx"],
                include_dirs=sys.path + ([os.path.dirname(gmpy2.__file__)] if ON_WINDOWS else []),
                library_dirs=sys.path + ([os.path.dirname(gmpy2.__file__)] if ON_WINDOWS else []),
                libraries=['mpc','mpfr','gmp']
            )
]

setup(
    name="cython_gmpy_test",
    ext_modules=cythonize(extensions, include_path=sys.path,
                          compiler_directives={'language_level' : "3"})
)
