from setuptools import Extension, setup
from Cython.Build import cythonize
import sys
import os
import platform
import gmpy2

gmpy2_packagedir = os.path.dirname(gmpy2.__file__)
library_dirs = sys.path + [gmpy2_packagedir]
libnames = ['mpc','mpfr','gmp']

bundled_libs = os.path.join(gmpy2_packagedir, '..', 'gmpy2.libs')
if os.path.isdir(bundled_libs):
    library_dirs += [bundled_libs]
    if platform.system() == 'Linux':
        libnames = [':' + d for d in os.listdir(bundled_libs)]
    elif platform.system() == 'Darwin':
        libnames = [':' + bundled_libs + d for d in os.listdir(bundled_libs)]


extensions = [
    Extension("test_cython", ["test_cython.pyx"],
              include_dirs=sys.path + [gmpy2_packagedir],
              library_dirs=library_dirs,
              libraries=libnames)]

setup(
    name="cython_gmpy_test",
    ext_modules=cythonize(extensions, include_path=sys.path,
                          compiler_directives={'language_level' : "3"})
)
