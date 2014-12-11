import sys
import os
from distutils.core import setup, Extension
from distutils.command.clean import clean
from distutils.command.build_ext import build_ext
from distutils.sysconfig import get_python_inc, get_python_lib

def writeln(s):
    sys.stdout.write('%s\n' % s)
    sys.stdout.flush()

# Some operating systems may use a different library directory under the
# prefix specified by --prefix. It must be manually changed.

lib_path = '/lib'

# Extract the version from MPFR/MPC

def get_mpfr_version(fname):
    result = []
    try:
        with open(fname) as f:
            for line in f:
                if line.startswith('#define MPFR_VERSION_MAJOR'):
                    result.append(int(line.split()[-1]))
                if line.startswith('#define MPFR_VERSION_MINOR'):
                    result.append(int(line.split()[-1]))
                if line.startswith('#define MPFR_VERSION_PATCHLEVEL'):
                    result.append(int(line.split()[-1]))
        return tuple(result)
    except:
        return (0,0,0)

def get_mpc_version(fname):
    result = []
    try:
        with open(fname) as f:
            for line in f:
                if line.startswith('#define MPC_VERSION_MAJOR'):
                    result.append(int(line.split()[-1]))
                if line.startswith('#define MPC_VERSION_MINOR'):
                    result.append(int(line.split()[-1]))
                if line.startswith('#define MPC_VERSION_PATCHLEVEL'):
                    result.append(int(line.split()[-1]))
        return tuple(result)
    except:
        return (0,0,0)

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

    # Extract the version information from the various header files. Since header
    # store the information differently, a separate function is provided for each
    # library.

    def check_versions(self):
        # Check the specified list of include directories to verify that valid
        # versions of MPFR and MPC are available. If so, add entries to the
        # appropriate lists

        # Find the directory specfied for PREFIX.
        prefix = ''
        for i,d in enumerate(self.extensions[0].define_macros[:]):
            if d[0] == 'PREFIX':
                prefix = d[1]
                try:
                    self.extensions[0].define_macros.remove(d)
                except ValueError:
                    pass

        if sys.version.find('MSC') == -1:
            windows = False
            base_dir = ['/usr']
            addin_dirs = ['/usr/local']
        else:
            windows = True
            base_dir = []
            addin_dirs = []

        if prefix:
            search_dirs =  base_dir + addin_dirs + [prefix]
        else:
            search_dirs = base_dir + addin_dirs

        if 'gmp' in self.extensions[0].libraries:
            mplib = 'gmp'
        else:
            mplib = 'mpir'

        use_mpfr = 'mpfr' in self.extensions[0].libraries
        use_mpc = 'mpc' in self.extensions[0].libraries

        # Try to find a directory prefix that contains valid MPFR and MPC
        # libraries. Only the version numbers of MPFR and MPC are checked.
        # Don't bother checking for GMP version since it is effectively a
        # prerequisite for MPFR and MPC.

        if not search_dirs:
            return

        mpfr_found = ''
        mpfr_version = (0,0,0)
        mpc_found = ''
        mpc_version = (0,0,0)

        for adir in search_dirs:
            lookin = adir + '/include'

            # If MPFR and MPC support is required, verify that the header files
            # exist in the same directory. If not, generate an error message.
            # If header isn't found, go to the next directory.

            # For debugging information, uncomment the following lines.
            # writeln('looking in: %s' % lookin)
            # writeln('mpfr.h found: %s' % os.path.isfile(lookin + '/mpfr.h'))
            # writeln('mpfr.h version %s' % repr(get_mpfr_version(lookin + '/mpfr.h')))
            # writeln('mpc.h found: %s' % os.path.isfile(lookin + '/mpc.h'))
            # writeln('mpc.h version %s' % repr(get_mpc_version(lookin + '/mpc.h')))

            if use_mpfr and os.path.isfile(lookin + '/mpfr.h'):
                v = get_mpfr_version(lookin + '/mpfr.h')
                if v >= mpfr_version:
                    mpfr_found = adir
                    mpfr_version = v

            if use_mpc and os.path.isfile(lookin + '/mpc.h'):
                v = get_mpc_version(lookin + '/mpc.h')
                if v >=mpc_version:
                    mpc_found = adir
                    mpc_version = v

        if mpfr_version < (3,1,0):
            mpfr_found = False
        if mpc_version < (1,0,0):
            mpc_found = False

        if (use_mpfr and not mpfr_found) or (use_mpc and not mpc_found):
            writeln('----------------------------------------------------------------')
            writeln('setup.py was not able to detect the required versions of MPFR')
            writeln('and/or MPC. gmpy2 requires MPFR version 3.1.0 or greater and')
            writeln('MPC version 1.0.0 or greater. To specify a directory prefix that')
            writeln('contains the proper versions, use the --prefix=<dir> option.')
            writeln('')
            writeln('In some circumstances, the correct versions may be present and')
            writeln('this warning can be ignored. If you have difficulties compiling')
            writeln('or running gmpy2, please try compiling with the --prefix option.')
            writeln('')
            writeln('It is possible to compile gmpy2 without support for MPFR and MPC')
            writeln('but that is not a supported configuration. Beginning with v2.1.0,')
            writeln('MPFR and MPC will be required.')
            writeln('')
            writeln('setup.py will continue and attempt to compile gmpy2.')
            writeln('-----------------------------------------------------------------')

        # Add the directory information for location where valid versions were
        # found. This can cause confusion if there are multiple installations of
        # the same version of Python on the system.

        for adir in (mpfr_found, mpc_found):
            if not adir:
                continue
            if adir in base_dir:
                continue
            if adir + '/include' in self.extensions[0].include_dirs:
                continue
            self.extensions[0].include_dirs += [adir + '/include']
            self.extensions[0].library_dirs += [adir + lib_path]
            if not windows:
                self.extensions[0].runtime_library_dirs += [adir + lib_path]

        # If the instance of Python used to compile gmpy2 not found in 'prefix',
        # then specify the Python shared library to use.
        #~ if not windows and not static and not get_python_lib(standard_lib=True).startswith(prefix):
           #~ self.extensions[0].libraries = [get_python_lib(standard_lib=True)] \
                                           #~ + self.extensions[0].libraries

        # For debugging information, uncomment the following lines.
        # writeln([mpfr_found, mpc_found])
        # writeln(self.extensions[0].include_dirs)
        # writeln(self.extensions[0].library_dirs)
        # writeln(self.extensions[0].runtime_library_dirs)
        # writeln(self.extensions[0].libraries)

    def finalize_options(self):
        build_ext.finalize_options(self)
        gmpy_build_ext.check_versions(self)
        # Check if --force was specified.
        for i,d in enumerate(self.extensions[0].define_macros[:]):
            if d[0] == 'FORCE':
                self.force = 1
                try:
                    self.extensions[0].define_macros.remove(d)
                except ValueError:
                    pass

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
#  --lib64         -> use /prefix/lib64 instead of /prefix/lib
#  --prefix=<...>  -> add the specified directory prefix to the beginning of
#                     the list of directories that are searched for GMP, MPFR,
#                     and MPC
#  --static        -> create a statically linked library
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

# If 'clean' is the only argument to setup.py then we want to skip looking for
# header files.

if sys.argv[1].lower() in ['build', 'install']:
    do_search = True
else:
    do_search = False

# Parse command line arguments. If custom prefix location is specified, it is
# passed as a define so it can be processed in the custom build_ext defined
# above.

defines = []

# Beginning with v2.1.0, MPFR and MPC will be required.

use_mpc = True
use_mpfr = True
force = False
static = False

for token in sys.argv[:]:
    if token.lower() == '--force':
        defines.append( ('FORCE', 1) )
        sys.argv.remove(token)

    if token.lower() == '--lib64':
        lib_path = '/lib64'
        sys.argv.remove(token)

    if token.lower() == '--mpir':
        mplib = 'mpir'
        sys.argv.remove(token)

    if token.lower() == '--gmp':
        mplib = 'gmp'
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
            defines.append( ('PREFIX', token.split('=')[1]) )

        except:
            writeln('Please include a directory location.')
        sys.argv.remove(token)

    if token.lower() == '--static':
        static = True
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
            defines.append( ('PREFIX', token.split('=')[1]) )
        except:
            writeln('Please include a directory location.')
        sys.argv.remove(token)
        writeln('The -DDIR option is deprecated. Use --prefix instead.')

incdirs = ['./src']
libdirs = []
rundirs = []
extras = []

# Specify extra link arguments for Windows.

if sys.version.find('MSC') == -1:
    my_extra_link_args = None
else:
    my_extra_link_args = ["/MANIFEST"]

mp_found = False

prefix = ''
for i,d in enumerate(defines):
    if d[0] == 'PREFIX':
        prefix = d[1]

if mplib == 'mpir':
    defines.append( ('MPIR', None) )
    libs = ['mpir']
    if static:
        extras.append(prefix + lib_path + '/libmpir.a')
else:
    libs = ['gmp']
    if static:
        extras.append(prefix + lib_path + '/libgmp.a')

if use_mpfr:
    defines.append( ('WITHMPFR', None) )
    libs.append('mpfr')
    if static:
        extras.append(prefix + lib_path + '/libmpfr.a')

if use_mpc:
    defines.append( ('WITHMPC', None) )
    libs.append('mpc')
    if static:
        extras.append(prefix + lib_path + '/libmpc.a')

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'

# prepare the extension for building

my_commands = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext}

gmpy2_ext = Extension('gmpy2',
                      sources=['src/gmpy2.c'],
                      include_dirs=incdirs,
                      library_dirs=libdirs,
                      libraries=libs,
                      runtime_library_dirs=rundirs,
                      define_macros = defines,
                      extra_objects = extras,
                      extra_link_args = my_extra_link_args)

setup(name = "gmpy2",
      version = "2.0.5",
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
