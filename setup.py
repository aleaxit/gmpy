import platform
import os
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

ON_WINDOWS = platform.system() == 'Windows'
_comp_args = ["DSHARED=1"]
sources = ['src/gmpy2.c']

# Utility function to read the contents of the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

class Gmpy2Build(build_ext):
    description = "Build gmpy2 with custom build options"
    user_options = build_ext.user_options + [
        ('fast', None,
         "Depend on MPFR and MPC internal implementations details"
         "(even more than the standard build)"),
        ('gcov', None, "Enable GCC code coverage collection"),
        ('vector', None, "Include the vector_XXX() functions;"
         "they are unstable and under active development"),
        ('static', None, "Enable static linking compile time options."),
        ('static-dir=', None, "Enable static linking and specify location."),
        ('gdb', None, "Build with debug symbols."),
    ]

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.fast = False
        self.gcov = False
        self.vector = False
        self.static = False
        self.static_dir = False
        self.gdb = False

    def finalize_options(self):
        build_ext.finalize_options(self)
        if self.fast:
            _comp_args.append('DFAST=1')
        if self.gcov:
            if ON_WINDOWS:
                raise ValueError("Cannot enable GCC code coverage on Windows")
            _comp_args.append('DGCOV=1')
            _comp_args.append('O0')
            _comp_args.append('-coverage')
            self.libraries.append('gcov')
        if self.vector:
            _comp_args.append('DVECTOR=1')
        if self.static:
            _comp_args.remove('DSHARED=1')
            _comp_args.append('DSTATIC=1')
        if self.gdb:
            _comp_args.append('ggdb')
        if self.static_dir:
            _comp_args.remove('DSHARED=1')
            _comp_args.append('DSTATIC=1')
            self.include_dirs.append(self.static_dir + '/include')
            self.library_dirs.append(self.static_dir + '/lib')

    def build_extensions(self):
        compiler = self.compiler.compiler_type
        _prefix = '-' if compiler != 'msvc' else '/'
        for i in range(len(_comp_args)):
            _comp_args[i] = ''.join([_prefix, _comp_args[i]])
        build_ext.build_extensions(self)

extensions = [
    Extension('gmpy2.gmpy2',
              sources=sources,
              include_dirs=['./src'],
              libraries=['mpc','mpfr','gmp'],
              extra_compile_args=_comp_args,
              )
]

cmdclass = {'build_ext': Gmpy2Build}

setup(
    name="gmpy2",
    version="2.1.0",
    author="Case Van Horsen",
    author_email="casevh@gmail.com",
    cmdclass=cmdclass,
    license="LGPL-3.0+",
    url="https://github.com/aleaxit/gmpy",
    description="gmpy2 interface to GMP/MPIR, MPFR, "
    "and MPC for Python 2.7 and 3.5+",
    long_description=read('README'),
    zip_safe=False,
    include_package_data=True,
    package_data={'gmpy2': [
        '*.pxd',
        'gmpy2.h',
        '*.dll',
    ]},
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: C',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    keywords="gmp mpfr mpc multiple-precision arbitrary-precision precision bignum",
    ext_modules=extensions,
)
