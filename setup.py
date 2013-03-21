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

# Improved clean command.

class gmpy_clean(clean):

    def run(self):
        self.all = True
        clean.run(self)

# Define a custom build class to force a new build.

class gmpy_build_ext(build_ext):

    def finalize_options(self):
        build_ext.finalize_options(self)
        self.force = 1

# Extract the version information from the various header files. Since header
# store the information differently, a separate function is provided for each
# library.

def get_mpfr_version(fname):
    result = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#define MPFR_VERSION_MAJOR'):
                result.append(int(line.split()[-1]))
            if line.startswith('#define MPFR_VERSION_MINOR'):
                result.append(int(line.split()[-1]))
            if line.startswith('#define MPFR_VERSION_PATCHLEVEL'):
                result.append(int(line.split()[-1]))
    return tuple(result)

def get_mpc_version(fname):
    result = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#define MPC_VERSION_MAJOR'):
                result.append(int(line.split()[-1]))
            if line.startswith('#define MPC_VERSION_MINOR'):
                result.append(int(line.split()[-1]))
            if line.startswith('#define MPC_VERSION_PATCHLEVEL'):
                result.append(int(line.split()[-1]))
    return tuple(result)

# Several command line options can be used to modify compilation of GMPY2. To
# maintain backwards compatibility with older versions of setup.py, the old
# options are still supported.
#
# New-style options
#
#  --force         -> ignore timestamps and recompile
#  --mpir          -> use MPIR instead of GMP (GMP is the default on
#                     non-Windows operating systems)
#  --gmp           -> use GMP instead of MPIR
#  --nompfr        -> disable MPFR and MPC library support
#  --nompc         -> disable MPC support (MPFR should still work)
#  --prefix=<...>  -> add the specified directory prefix to the beginning of
#                     the list of directories that are searched for GMP, MPFR,
#                     and MPC
#
# Old-stype options
#
#   -DMPIR       -> use MPIR instead of GMP
#   -DGMP        -> use GMP instead of MPIR
#   -DNOMPFR     -> disable MPFR and MPC library support
#   -DNOMPC      -> disable MPC support (MPFR should still work)
#   -Ddir=<...>  -> add the specified directory to beginning of the list of
#                   directories that are searched for GMP, MPFR, and MPC

# Windows build defaults to using MPIR.

if sys.version.find('MSC') == -1:
    mplib='gmp'
else:
    mplib='mpir'

# Specify the default search directories for Unix/Linux/MacOSX. The prefixes
# '/usr/local' and '/usr' are searched since they are "standard" on most Linux
# distributions.
#
# The directories are searched for 'include/gmp.h' or 'include/mpir.h'. If
# found, then 'include/mpfr.h' and 'include/mpc.h' are expected to be found
# under the same prefix.
#
# To specify an alternate location for the gmp/mpir, mpfr, and mpc, use the
# --prefix=<<prefix>> option. If --prefix is specified, the standard locations
# are not searched. When using --prefix, the run_path option is specified when
# linking gmpy2. The is done to make it easy to support locally compiled
# versions of GMP/MPIR, MPFR, and MPC.

if sys.version.find('MSC') == -1:
    search_dirs = ['/usr/local', '/usr']
else:
    search_dirs = []

# Some operating systems may use a different library directory under the
# prefix specified by --prefix. It must be manually changed.

lib_path = '/lib'

# If 'clean' is the only argument to setup.py then we want to skip looking for
# header files.

if len(sys.argv) == 2 and sys.argv[1].lower() == 'clean':
    do_search = False
else:
    do_search = True

use_mpc = True
use_mpfr = True
force = False
prefix = False

for token in sys.argv[:]:
    if token.lower() == '--force':
        force = True
        sys.argv.remove(token)

    if token.lower() == '--mpir':
        mplib='mpir'
        sys.argv.remove(token)

    if token.lower() == '--gmp':
        mplib='gmp'
        sys.argv.remove(token)

    if token.lower() == '--nompc':
        use_mpc = False
        sys.argv.remove(token)

    if token.lower() == '--nompfr':
        use_mpfr = False
        use_mpc = False
        sys.argv.remove(token)

    if token.lower().startswith('--prefix'):
        try:
            prefix = True
            search_dirs = [token.split('=')[1]]
        except:
            writeln('Please include a directory location.')
        sys.argv.remove(token)

    # The following options are deprecated and will be removed in the future.
    if token.upper().startswith('-DMPIR'):
        mplib='mpir'
        sys.argv.remove(token)
        writeln('The -DMPIR option is deprecated. Use --mpir instead.')

    if token.upper().startswith('-DGMP'):
        mplib='gmp'
        sys.argv.remove(token)
        writeln('The -DGMP option is deprecated. Use --gmp instead.')

    if token.upper().startswith('-DNOMPC'):
        use_mpc = False
        sys.argv.remove(token)
        writeln('The -DNOMPC option is deprecated. Use --nompc instead.')

    if token.upper().startswith('-DNOMPFR'):
        use_mpfr = False
        use_mpc = False
        sys.argv.remove(token)
        writeln('The -DNOMPFR option is deprecated. Use --nompfr instead.')

    if token.upper().startswith('-DDIR'):
        try:
            prefix = True
            search_dirs = [token.split('=')[1]]
        except:
            writeln('Please include a directory location.')
        sys.argv.remove(token)
        writeln('The -DDIR option is deprecated. Use --prefix instead.')

incdirs = ['./src']
libdirs = []
rundirs = []

# Specify extra link arguments for Windows.

if sys.version.find('MSC') == -1:
    my_extra_link_args = None
else:
    my_extra_link_args = ["/MANIFEST"]

mp_found = False

# TODO: Parse the header files and check for the appropriate version.

if do_search:
    if not search_dirs:
        writeln('Please specify the prefix directory for the include and library files')
        writeln('using the --prefix=<dir> option.');
        sys.exit()

    for adir in search_dirs:
        lookin = adir + '/include'
        if os.path.isfile(lookin + '/' + mplib + '.h'):
            mp_found = True
            # Only modify the inc/lib/run directories if --prefix was used.
            if prefix:
                incdirs += [lookin]
                libdirs += [adir + lib_path]
                rundirs = [adir + lib_path]

            # If MPFR and MPC support is required, verify that the header files
            # exist in the same directory. If not, generate an error message.
            if use_mpfr and not os.path.isfile(lookin + '/mpfr.h'):
                writeln('mpfr.h is not present in %s.' % lookin)
                writeln('To disable support for MPFR, use the --nompfr option.')
                writeln('To specify a directory prefix for the include and library files,')
                writeln('use the --prefix=<dir> option.')
                sys.exit()
            if use_mpc and not os.path.isfile(lookin + '/mpc.h'):
                writeln('mpc.h is not present in %s.' % lookin)
                writeln('To disable support for MPC, use the --nompc option.')
                writeln('To specify a directory prefix for the include and library files,')
                writeln('use the --prefix=<dir> option.')
                sys.exit()

            # Check for the proper versions of MPFR and MPC.
            mpfr_version = get_mpfr_version(lookin + '/mpfr.h')
            if use_mpfr and mpfr_version < (3,1,0):
                writeln('MPFR version %s.%s.%s was found.' % mpfr_version)
                writeln('The mininum required version is 3.1.0.')
                writeln('To specify the location of an updated version, use the')
                writeln('--prefix=<dir> option.')
                sys.exit()
            mpc_version = get_mpc_version(lookin + '/mpc.h')
            if use_mpc and mpc_version < (1,0,0):
                writeln('MPC version %s.%s.%s was found.' % mpc_version)
                writeln('The mininum required version is 1.0.0.')
                writeln('To specify the location of an updated version, use the')
                writeln('--prefix=<dir> option.')
                sys.exit()

            break

    if not mp_found:
        writeln('%s.h could not be found.' % mplib)
        writeln('To specify a directory prefix for the include and library files,')
        writeln('use the --prefix=<dir> option.')
        sys.exit()

# Configure the defines...

defines = []
if mplib == 'mpir':
    defines.append( ('MPIR', 1) )
if use_mpfr:
    defines.append( ('WITHMPFR', 1) )
if use_mpc:
    defines.append( ('WITHMPC', 1) )

# Build list of the required libraries. If the instance of Python used to compile
# gmpy2 not found under --prefix, then specify the Python shared to use.

if prefix:
    if get_python_lib(standard_lib=True).startswith(search_dirs[0]):
        libs = [mplib]
    else:
        libs = [get_python_lib(standard_lib=True), mplib]
else:
    libs = [mplib]

if use_mpfr:
    libs.append('mpfr')
if use_mpc:
    libs.append('mpc')

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'

# prepare the extension for building

if force:
    cmdclass = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext}
else:
    cmdclass = {'clean' : gmpy_clean}

gmpy2_ext = Extension('gmpy2',
                      sources=['src/gmpy2.c'],
                      include_dirs=incdirs,
                      library_dirs=libdirs,
                      libraries=libs,
                      runtime_library_dirs=rundirs,
                      define_macros = defines,
                      extra_link_args = my_extra_link_args)

setup(name = "gmpy2",
      version = "2.0.0b5",
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
      cmdclass = cmdclass,
      ext_modules = [gmpy2_ext]
)
