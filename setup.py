import platform
import os
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.command.install_data import install_data

ON_WINDOWS = platform.system() == 'Windows'
_comp_args = ["DSHARED=1"]
link_args = []

sources = ['src/gmpy2.c']
_libs = ['mpfr', 'mpc']

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
        ('mpir', None, "Enable use of mpir library instead of gmp."
         "gmp is the default on Posix systems while mpir the default on"
         "Windows and MSVC"),
        ('static', None, "Enable static linking compile time options."),
        ('gdb', None, "Build with debug symbols."),
    ]

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.fast = False
        self.gcov = False
        self.vector = False
        self.mpir = False
        self.static = False
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
            _comp_args.append('--coverage')
            link_args.append('--coverage')
            _libs.append('gcov')
        if self.vector:
            _comp_args.append('DVECTOR=1')
        if self.static:
            _comp_args.remove('DSHARED=1')
            _comp_args.append('DSTATIC=1')
        if self.gdb:
            _comp_args.append('ggdb')

    def build_extensions(self):
        compiler = self.compiler.compiler_type
        if compiler == 'mingw32':
            _libs.append('mpfr')
            _libs.append('mpc')
            _comp_args.append('DMSYS2=1')
            if self.mpir:
                _comp_args.append('DMPIR=1')
                _libs.append('mpir')
            else:
                _libs.append('gmp')
        elif self.mpir or ON_WINDOWS:
            # --mpir or on Windows and MSVC
            _comp_args.append('DMPIR=1')
            _libs.append('mpir')
            if ON_WINDOWS and not self.static:
                # MSVC shared build
                _comp_args.append('MSC_USE_DLL')
        else:
            _libs.append('gmp')
        _prefix = '-' if compiler != 'msvc' else '/'
        for i in range(len(_comp_args)):
            _comp_args[i] = ''.join([_prefix, _comp_args[i]])
        build_ext.build_extensions(self)


extensions = [
    Extension('gmpy2.gmpy2',
              sources=sources,
              include_dirs=['./src'],
              libraries=_libs,
              extra_compile_args=_comp_args,
              extra_link_args=link_args,
              )
]

class gmpy_install_data(install_data):
    def finalize_options(self):
        install_data.finalize_options(self)
        install = self.distribution.get_command_obj('install')
        self.install_dir = install.install_purelib

cmdclass = {'build_ext': Gmpy2Build, 'install_data' : gmpy_install_data}

setup(
    name="gmpy2",
    version="2.1.0a3dev0",
    author="Case Van Horsen",
    author_email="casevh@gmail.com",
    cmdclass=cmdclass,
    license="LGPL-3.0+",
    url="https://github.com/aleaxit/gmpy",
    description="gmpy2 interface to GMP/MPIR, MPFR, "
    "and MPC for Python 2.6+ and 3.4+",
    long_description=read('README'),
    zip_safe=False,
    data_files = [('', ['src/gmpy2.pxd']), ('', ['src/gmpy2.h'])],
    packages=find_packages(),
    classifiers=[
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
