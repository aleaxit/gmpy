# pysymbolicext is no longer supported. It is included soley to test the
# exported C API.
import sys
from distutils.core import setup, Extension

if sys.version.find('MSC')==-1:
    gmpy_ext = Extension('pysymbolicext', sources=['src/pysymbolicext.c'],
        libraries=['gmp'])
else:
    gmpy_ext = Extension('pysymbolicext', sources=['src/pysymbolicext.c'],
        libraries=['gmp'],include_dirs=['./src'])

setup (name = "pysymbolicext",
       version = "0.2",
       description = "PySymbolic Python/GMP extensions (Pollard's rho)",
       author = "Pearu Peterson",
       maintainer = "Alex Martelli",
       maintainer_email = "aleaxit@gmail.com",
       url = "http://code.google.com/p/gmpy/source/",

       ext_modules = [ gmpy_ext ]
)
