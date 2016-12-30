# gmpy2_demo is not supported. It is included soley to test the
# exported C API.

import sys
from distutils.core import setup, Extension

if sys.version.find('MSC')==-1:
    gmpy_ext = Extension('gmpy2_demo', sources=['gmpy2_demo.c'],
        libraries=['gmp'])
else:
    gmpy_ext = Extension('gmpy2_demo', sources=['gmpy2_demo.c'],
        libraries=['gmp'],include_dirs=['.'])

setup (name = "gmpy2_demo",
       version = "0.3",
       description = "gmpy2_demo: gmpy2 demonstration programs",
       author = "Case Van Horsen",
       maintainer = "Case Van Horsen",
       maintainer_email = "casevh@gmail.com",
       url = "https://github.com/aleaxit/gmpy",

       ext_modules = [ gmpy_ext ]
)
