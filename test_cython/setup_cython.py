from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys

extensions = [
    Extension("test_cython", ["test_cython.pyx"],
        include_dirs=sys.path,
        )
]

setup(
    name="cython_gmpy_test",
    ext_modules=cythonize(extensions, include_path=sys.path,
                          compiler_directives={'language_level' : "3"})
)
