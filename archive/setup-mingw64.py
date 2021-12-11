# This file is intended soley to build statically-linked binaries for Windows.
#
# Please read the file "mingw64_build.txt" for information for required
# environment.
#
import sys

if sys.platform != "win32":
    raise ValueError("This setup file only runs on Windows.")

def writeln(s):
    sys.stdout.write('%s\n' % s)
    sys.stdout.flush()

# Utility function to read the contents of the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# Fail gracefully for old versions of Python.
# Note: The slice will return an incorrect value beginning with Python 3.10
#       but the result of the comparison is still valid.

if sys.version[:3] < '2.6':
    writeln("GMPY2 requires Python 2.7 or later.")
    writeln("Please use GMPY 1.x for earlier versions of Python.")
    sys.exit()

from distutils.core import setup, Extension
from distutils.command.clean import clean
from distutils.command.build_ext import build_ext
from distutils.command.install_data import install_data

# Run python setup.py build_ext --help for build options

# Improved clean command.

class gmpy_clean(clean):

    def run(self):
        self.all = True
        clean.run(self)

build_options = [
    ('fast',    None, 'depend on MPFR and MPC internal implementations details'),
    ('vector',  None, 'include the vector_XXX() functions; they are unstable and under active development'),
    ('shared',  None, 'Build using shared libraries'),
    ('static',  None, 'Build using static libraries'),
]

# Define a custom build class that parses to the defined macros to alter the
# build setup.

class gmpy_build_ext(build_ext):
    default_opts = build_ext.user_options
    user_options = list(build_options)
    user_options.extend(default_opts)

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.fast = False
        self.vector = False
        self.static = None
        self.shared = None

    def doit(self):

        defines = []

        if self.static is None and self.shared is None:
            self.shared = True

        if self.vector:
            defines.append(('VECTOR', 1))

        if self.fast:
            defines.append(('FAST', 1))

        if self.shared:
            defines.append(('SHARED', 1))

        if self.static:
            defines.append(('STATIC', 1))

        self.extensions[0].libraries.extend(['gmp', 'mpfr', 'mpc'])

        self.extensions[0].define_macros.extend(defines)

    def finalize_options(self):
        build_ext.finalize_options(self)
        gmpy_build_ext.doit(self)

# custom install data in order that data_files
# get installed together with the .so

class gmpy_install_data(install_data):
    def finalize_options(self):
        install_data.finalize_options(self)
        install = self.distribution.get_command_obj('install')
        self.install_dir = install.install_purelib

# prepare the extension for building

my_commands = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext, 'install_data' : gmpy_install_data}

gmpy2_ext = Extension('gmpy2',
                      sources=[os.path.join('src', 'gmpy2.c')],
                      include_dirs=['./src'])

setup(name = "gmpy2",
      version = "2.1.0b5",
      author = "Case Van Horsen",
      author_email = "casevh@gmail.com",
      license = "LGPL-3.0+",
      url = "https://github.com/aleaxit/gmpy",
      description = "gmpy2 interface to GMP, MPFR, and MPC for Python 2.7 and 3.5+",
      long_description = read('README'),
      data_files = [('', ['src/gmpy2.pxd']), ('gmpy2', ['src/gmpy2.h'])],
      classifiers = [
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
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      keywords = "gmp mpir mpfr mpc multiple-precision arbitrary-precision precision bignum",
      cmdclass = my_commands,
      ext_modules = [gmpy2_ext]
)
