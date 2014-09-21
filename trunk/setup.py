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
# prefix specified by --prefix. It can be changed using the --lib64 option.

lib_path = 'lib'

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
            search_dirs = base_dir + addin_dirs + [prefix]
        else:
            search_dirs = base_dir + addin_dirs

        if 'gmp' in self.extensions[0].libraries:
            mplib = 'gmp'
        else:
            mplib = 'mpir'

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
            lookin = os.path.join(adir, 'include')

            # For debugging information, uncomment the following lines.
            # writeln('looking in: %s' % lookin)
            # writeln('mpfr.h found: %s' % os.path.isfile(os.path.join(lookin, 'mpfr.h')))
            # writeln('mpfr.h version %s' % repr(get_mpfr_version(os.path.join(lookin, 'mpfr.h'))))
            # writeln('mpc.h found: %s' % os.path.isfile(os.path.join(lookin, 'mpc.h')))
            # writeln('mpc.h version %s' % repr(get_mpc_version(os.path.join(lookin, 'mpc.h'))))

            if os.path.isfile(os.path.join(lookin, 'mpfr.h')):
                v = get_mpfr_version(os.path.join(lookin, 'mpfr.h'))
                if v >= mpfr_version:
                    mpfr_found = adir
                    mpfr_version = v

            if os.path.isfile(os.path.join(lookin, 'mpc.h')):
                v = get_mpc_version(os.path.join(lookin, 'mpc.h'))
                if v >= mpc_version:
                    mpc_found = adir
                    mpc_version = v

        if mpfr_version < (3,1,0):
            mpfr_found = False
        if mpc_version < (1,0,0):
            mpc_found = False

        if not mpfr_found or not mpc_found:
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
            writeln('setup.py will still continue and attempt to compile gmpy2.')
            writeln('-----------------------------------------------------------------')

        # Add the directory information for location where valid versions were
        # found. This can cause confusion if there are multiple installations of
        # the same version of Python on the system.

        for adir in (mpfr_found, mpc_found):
            if not adir:
                continue
            if adir in base_dir:
                continue
            if os.path.join(adir, 'include') in self.extensions[0].include_dirs:
                continue
            self.extensions[0].include_dirs += [os.path.join(adir, 'include')]
            self.extensions[0].library_dirs += [os.path.join(adir, lib_path)]
            if not windows and adir not in base_dir:
                self.extensions[0].runtime_library_dirs += [os.path.join(adir, lib_path)]

        # If the instance of Python used to compile gmpy2 not found in 'prefix',
        # then specify the Python shared library to use.
        #~ if not windows and not get_python_lib(standard_lib=True).startswith(prefix):
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
        # Check if --msys2 was specified.
        for i,d in enumerate(self.extensions[0].define_macros[:]):
            if d[0] == 'MSYS2':
                self.compiler = 'mingw32'

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
#  --msys2         -> build on Windows using MSYS2, MinGW, and GMP
#  --nompfr        -> was disable MPFR and MPC library support in 2.0.x
#                     MPFR and MPC are required for 2.1.x
#  --nompc         -> was disable MPC support in 2.0.x
#                     MPFR and MPC are required 2.1.x
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
#   -DNOMPFR     -> was disable MPFR and MPC library support in 2.0.x
#                   MPFR and MPC are required for 2.1.x
#   -DNOMPC      -> was disable MPC support (MPFR should still work)
#                   MPFR and MPC are required for 2.1.x
#   -Ddir=<...>  -> add the specified directory to beginning of the list of
#                   directories that are searched for GMP, MPFR, and MPC

# Windows build defaults to using MPIR. This will be over-written later if
# --msys2 is given.

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

force = False
static = False
msys2 = False

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

    if token.lower() == '--msys2':
        msys2 = True
        mplib = 'gmp'
        static = True
        defines.append( ('MSYS2', 1) )
        sys.argv.remove(token)

    if token.lower() == '--nompc':
        writeln('--nompc is no longer supported. MPC is required.')
        sys.argv.remove(token)

    if token.lower() == '--nompfr':
        writeln('--nompfr is no longer supported. MPFR is required.')
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
        writeln('-DNOMPC is no longer supported. MPC is required.')
        sys.argv.remove(token)

    if token.upper().startswith('-DNOMPFR'):
        writeln('-DNOMPFR is no longer supported. MPFR is required.')
        sys.argv.remove(token)

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

if sys.version.find('MSC') == -1 or msys2:
    my_extra_link_args = None
else:
    my_extra_link_args = ["/MANIFEST"]

mp_found = False

# Configure the defines...

prefix = ''
for i,d in enumerate(defines):
    if d[0] == 'PREFIX':
        prefix = d[1]

if mplib == 'mpir':
    defines.append( ('MPIR', None) )
    libs = ['mpir']
    if static:
        extras.append(os.path.join(prefix, lib_path, 'libmpir.a'))
else:
    libs = ['gmp']
    if static:
        extras.append(os.path.join(prefix, lib_path, 'libgmp.a'))

libs.append('mpfr')
if static:
    extras.append(os.path.join(prefix, lib_path, 'libmpfr.a'))

libs.append('mpc')
if static:
    extras.append(os.path.join(prefix, lib_path, 'libmpc.a'))

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'

# prepare the extension for building

my_commands = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext}

gmpy2_ext = Extension('gmpy2',
                      sources=[os.path.join('src', 'gmpy2.c')],
                      include_dirs=incdirs,
                      library_dirs=libdirs,
                      libraries=libs,
                      runtime_library_dirs=rundirs,
                      define_macros = defines,
                      extra_objects = extras,
                      extra_link_args = my_extra_link_args)

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
