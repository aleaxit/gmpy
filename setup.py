import platform
from setuptools import setup, find_packages, Extension

cpython = platform.python_implementation() == 'CPython'

sources = ['src/gmpy2.c']
_libs = ['gmp', 'mpfr', 'mpc'] if platform.system() != 'Windows' else [
    'libmpir', 'libmpfr', 'libmpc']
#     'Ws2_32', 'user32']

# _comp_args = ["-ggdb"]
_comp_args = ["-DSHARED=1"] if platform.system() != 'Windows' else [
    '-DMPIR=1']
extensions = [
    Extension('gmpy2/_gmpy2',
              sources=sources,
              include_dirs=["/usr/local/include", './src'],
              libraries=_libs,
              extra_compile_args=_comp_args,
              # extra_link_args=<..>,
          )]

setup(name="gmpy2",
      version="2.1.0a1",
      author="Case Van Horsen",
      author_email="casevh@gmail.com",
      license="LGPL-3.0+",
      url="https://github.com/aleaxit/gmpy",
      description="gmpy2 interface to GMP/MPIR, MPFR, and MPC for Python 2.6+ and 3.4+",
      include_package_data=True,
      package_data={'': [
          'src/gmpy2.pxd',
          'src/gmpy2.h',
      ]},
      packages=find_packages(),
      classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: C',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      keywords="gmp mpir mpfr mpc multiple-precision arbitrary-precision precision bignum",
      ext_modules=extensions,
)
