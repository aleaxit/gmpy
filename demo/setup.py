# gmpy2_demo is not supported. It is included soley to test the
# exported C API.

from setuptools import Extension, setup
import sys
import os
import platform
import gmpy2

gmpy2_packagedir = os.path.dirname(gmpy2.__file__)
library_dirs = sys.path + [gmpy2_packagedir]
libnames = ['mpc','mpfr','gmp']

if platform.system() != 'Windows':
    bundled_libs = gmpy2_packagedir+'/../gmpy2.libs/'
    if os.path.isdir(bundled_libs):
        library_dirs += [bundled_libs]
        if platform.system() == 'Linux':
            libnames = [':' + d for d in os.listdir(bundled_libs)]
        else:
            libnames = [':' + bundled_libs + d for d in os.listdir(bundled_libs)]

gmpy_ext = [
    Extension("gmpy2_demo", sources=["gmpy2_demo.c"],
              include_dirs=sys.path + [gmpy2_packagedir],
              library_dirs=library_dirs,
              libraries=libnames)]

setup (name = "gmpy2_demo",
       version = "0.3",
       description = "gmpy2_demo: gmpy2 demonstration programs",
       author = "Case Van Horsen",
       maintainer = "Case Van Horsen",
       maintainer_email = "casevh@gmail.com",
       url = "https://github.com/aleaxit/gmpy",

       ext_modules = gmpy_ext
)
