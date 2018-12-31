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
# prefix specified by --shared. It must be manually changed.

lib_path = '/lib'

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

        # Find the directory specfied for SHARED or STATIC.
        prefix = []
        for d in self.extensions[0].define_macros:
            if d[0] in ('SHARED', 'STATIC'):
                if d[1]:
                    prefix.extend(map(os.path.expanduser, d[1].split(":")))
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
            search_dirs = base_dir + addin_dirs + prefix
        else:
            search_dirs = base_dir + addin_dirs

        if 'gmp' in self.extensions[0].libraries:
            mplib = 'gmp'
        else:
            mplib = 'mpir'

        use_mpfr = 'mpfr' in self.extensions[0].libraries
        use_mpc = 'mpc' in self.extensions[0].libraries

        if not search_dirs:
            return

        gmp_found = ''
        mpfr_found = ''
        mpc_found = ''

        for adir in search_dirs:
            lookin = os.path.join(adir, 'include')

            if os.path.isfile(os.path.join(lookin, mplib + '.h')):
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
            if adir in base_dir:
                continue
            if os.path.join(adir, 'include') in self.extensions[0].include_dirs:
                continue
            self.extensions[0].include_dirs += [os.path.join(adir, 'include')]
            self.extensions[0].library_dirs += [os.path.join(adir, lib_path)]
            if shared and not windows:
                self.extensions[0].runtime_library_dirs += [os.path.join(adir, lib_path)]

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
#  --lib64         -> use /<...>/lib64 instead of /<...>/lib
#  --shared=<...>  -> add the specified directory prefix to the beginning of
#                     the list of directories that are searched for GMP, MPFR,
#                     and MPC shared libraries
#  --static=<...>  -> create a statically linked library using static files from
#                     the specified directory


# Windows build defaults to using MPIR.

if sys.version.find('MSC') == -1:
    mplib='gmp'
else:
    mplib='mpir'

# If 'clean' is the only argument to setup.py then we want to skip looking for
# header files.

if sys.argv[1].lower() in ['build', 'build_ext', 'install']:
    do_search = True
else:
    do_search = False

# Parse command line arguments. If custom prefix location is specified, it is
# passed as a define so it can be processed in the custom build_ext defined
# above.

defines = []

# Beginning with v2.1.0, MPFR and MPC will be required.

force = False
static = False
shared = False

for token in sys.argv[:]:
    if token.lower() == '--force':
        force = True
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

    if token.lower().startswith('--shared'):
        shared = True
        try:
            defines.append(('SHARED', token.split('=')[1]))
        except IndexError:
            defines.append(('SHARED', None))
        sys.argv.remove(token)

    if token.lower().startswith('--static'):
        static = True
        try:
            defines.append(('STATIC', token.split('=')[1]))
        except IndexError:
            defines.append(('STATIC', None))
        sys.argv.remove(token)

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

for i,d in enumerate(defines[:]):
    if d[0] in ('SHARED', 'STATIC'):
        if d[1]:
            prefix = d[1]
        defines.append((d[0], None))

writeln(prefix)

if force:
    defines.append(('FORCE', None))

if mplib == 'mpir':
    defines.append(('MPIR', None))
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

writeln(str(defines))

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'

# prepare the extension for building

my_commands = {'clean' : gmpy_clean, 'build_ext' : gmpy_build_ext}

gmpy2_ext = Extension('gmpy2',
                      sources=[os.path.join('src', 'gmpy2.c')],
                      libraries=libs,
                      define_macros = defines,
                      extra_objects = extras,
                      extra_link_args = my_extra_link_args)

setup(name = "gmpy2",
      version = "2.1.0a3",
      maintainer = "Case Van Horsen",
      maintainer_email = "casevh@gmail.com",
      url = "http://code.google.com/p/gmpy/",
      description = "GMP/MPIR, MPFR, and MPC interface to Python 2.6+ and 3.x",
      classifiers = [
        'Development Status :: 3 - Alpha',
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
