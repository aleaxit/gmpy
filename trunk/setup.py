import sys
import os
from distutils.core import setup, Extension
from distutils.command.clean import clean
from distutils.command.build_ext import build_ext
from distutils.sysconfig import get_python_inc, get_python_lib

def writeln(s):
    sys.stdout.write('%s\n' % s)
    sys.stdout.flush()


# Fail gracefully for old versions of Python.

if sys.version[:3] < '2.6':
    writeln("GMPY2 requires Python 2.6 or later.")
    writeln("Please use GMPY 1.x for earlier versions of Python.")
    sys.exit()

# Initialize some global values.

lib_path = 'lib'

if sys.version.find('MSC') == -1:
    windows = False
else:
    windows = True

# Several command line options can be used to modify compilation of GMPY2. 
#
#  --msys2         -> build on Windows using MSYS2, MinGW, and GMP
#  --lib64         -> use /prefix/lib64 instead of /prefix/lib
#  --lib32         -> use /prefix/lib32 instead of /prefix/lib
#  --shared=<...>  -> add the specified directory prefix to the beginning of
#                     the list of directories that are searched for GMP, MPFR,
#                     and MPC shared libraries
#  --static=<...>  -> create a statically linked library using libraries from
#                     specified path
#
# Ugly hack ahead. Sorry.
#
# I haven't found any examples on how to extend distutils with user-defined
# options. And the documentation is not helpful, either. So instead of fighting
# distutils, soem of the options are converted into macro definitions. Macros
# then get parsed by distutils and then a custom class reads the macros and
# tweaks the setup.
#
# The custom command line arguments are processed as follows:
#
#  --force is temporarily added to the list of defines and then removed by
#  gmpy_build_ext. It is used to allow setup.py install --force to work as
#  expected.
#
#  --msys2 is converted to -DMSYS2. Since MSYS2 needs to be defined as a macro
#  to control options in the GMPY2 source code, this make sense.
#
#  --lib64 and --lib32 are removed from sys.argv and the global variable
#  lib_path is set to 'lib64' or 'lib32' as appropriate.
#
#  --shared and --static are converted to -DSHARED and -DSTATIC.

defines = []

for token in sys.argv[:]:
    if token.lower() == '--force':
        defines.append( ('FORCE', 1) )
        sys.argv.remove(token)

    if token.lower() == '--lib64':
        lib_path = 'lib64'
        sys.argv.remove(token)

    if token.lower() == '--lib32':
        lib_path = 'lib32'
        sys.argv.remove(token)

    if token.lower() == '--msys2':
        defines.append( ('MSYS2', 1) )
        sys.argv.remove(token)

    if token.lower().startswith('--shared'):
        try:
            defines.append( ('SHARED', token.split('=')[1]) )
        except IndexError:
            pass
        sys.argv.remove(token)

    if token.lower().startswith('--static'):
        try:
            defines.append( ('STATIC', token.split('=')[1]) )
        except IndexError:
            pass
        sys.argv.remove(token)

# Improved clean command.

class gmpy_clean(clean):

    def run(self):
        self.all = True
        clean.run(self)

# Define a custom build class that parses to the defined macros to alter the
# build setup.

class gmpy_build_ext(build_ext):

    def initialize_options(self):
        build_ext.initialize_options(self)
        
    def doit(self):
        # Find the directory specfied for non-standard library location.
        search_dirs = []
        static = False
        msys2 = False

        # Assume that we will always want to use the GMP, MPFR, and MPC libraries.
        self.extensions[0].libraries.extend(['gmp', 'mpfr', 'mpc'])
        
        for d in self.extensions[0].define_macros[:]:
            if d[0] == 'MSYS2':
                self.compiler = 'mingw32'
                msys2 = True

            if d[0] == 'FORCE':
                self.force = 1
                try:
                    self.extensions[0].define_macros.remove(d)
                except ValueError:
                    pass
                    
            if d[0] in ('SHARED', 'STATIC'):
                if d[0] == 'STATIC':
                    static = True
                if d[0] == 'SHARED':
                    static = False
                if d[1]:
                    search_dirs.extend(map(os.path.expanduser, d[1].split(":")))
                    try:
                        self.extensions[0].define_macros.remove(d)
                    except ValueError:
                        pass

        # If non-default directories have been specified, we need to find the
        # exact location of the libraries to allow static or runtime linking.
        
        gmp_found = ''
        mpfr_found = ''
        mpc_found = ''
        if search_dirs:

            for adir in search_dirs:
                lookin = os.path.join(adir, 'include')
                if os.path.isfile(os.path.join(lookin, 'gmp.h')):
                    gmp_found = adir
                if os.path.isfile(os.path.join(lookin, 'mpfr.h')):
                    mpfr_found = adir
                if os.path.isfile(os.path.join(lookin, 'mpc.h')):
                    mpc_found = adir

        # Add the directory information for location where valid versions were
        # found. This can cause confusion if there are multiple installations of
        # the same version of Python on the system.

        for adir in (gmp_found, mpfr_found, mpc_found):
            if not adir:
                continue
            if os.path.join(adir, 'include') in self.extensions[0].include_dirs:
                continue
            self.extensions[0].include_dirs += [os.path.join(adir, 'include')]
            self.extensions[0].library_dirs += [os.path.join(adir, lib_path)]

            # Add the runtime linking options.
            if not static and not windows:
                self.extensions[0].runtime_library_dirs += [os.path.join(adir, lib_path)]

        # Add the static linking options.
        if static and gmp_found:
            self.extensions[0].extra_objects.append(os.path.join(gmp_found, lib_path, 'libgmp.a'))
        if static and mpfr_found:
            self.extensions[0].extra_objects.append(os.path.join(mpfr_found, lib_path, 'libmpfr.a'))
        if static and mpc_found:
            self.extensions[0].extra_objects.append(os.path.join(mpc_found, lib_path, 'libmpc.a'))

        # Add MSVC specific options.
        if windows and not msys2:
            self.extensions[0].extra_link_args.append('/MANIFEST')

    def finalize_options(self):
        build_ext.finalize_options(self)
        gmpy_build_ext.doit(self)
        

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'

# prepare the extension for building

my_commands = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext}

gmpy2_ext = Extension('gmpy2',
                      sources=[os.path.join('src', 'gmpy2.c')],
                      include_dirs=['./src'],
                      define_macros = defines)

setup(name = "gmpy2",
      version = "2.1.0a0",
      maintainer = "Case Van Horsen",
      maintainer_email = "casevh@gmail.com",
      url = "http://code.google.com/p/gmpy/",
      description = "GMP/MPIR, MPFR, and MPC interface to Python 2.6+ and 3.x",
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research'
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: C',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      cmdclass = my_commands,
      ext_modules = [gmpy2_ext]
)
