import sys, os
from distutils.core import setup, Extension

# Fail gracefully for old versions of Python.
if sys.version[:3] < '2.6':
    sys.stdout.write("GMPY2 requires Python 2.6 or later.\n")
    sys.stdout.write("Please use GMPY 1.x for earlier versions of Python.\n")
    sys.exit()

# Check for build options:
#   -DMPIR  -> use MPIR instead of GMP
#   -DGMP   -> use GMP instead of MPIR
# Windows build defaults to MPIR.
if sys.version.find('MSC') == -1:
    mplib='gmp'
else:
    mplib='mpir'

local_dir = None
running_build = False
for token in sys.argv:
    if token.upper().startswith('-DMPIR'):
        mplib='mpir'
    if token.upper().startswith('-DGMP'):
        mplib='gmp'
    if token.upper().startswith('BUILD'):
        running_build = True
    if token.upper().startswith('-DDIR'):
        try:
            local_dir = token.split('=')[1]
        except:
            pass

use_mpc = True
use_mpfr = True

incdirs = None
libdirs = None
rundirs = None
my_extra_link_args = None

# determine include and library dirs
if running_build and sys.version.find('MSC') == -1:
    # Unix-like build (including MacOSX)
    incdirs = ['./src']
    dirord = ['/opt/local', '/opt', '/usr/local']
    if local_dir:
        dirord = [local_dir] + dirord
    for adir in dirord:
        lookin = '%s/include' % adir
        if os.path.isfile(lookin + '/' + mplib + '.h'):
            incdirs = [lookin]
            # Verify that MPFR and MPC exist in the same directory
            if not os.path.isfile(lookin + '/mpfr.h'):
                use_mpfr = False
            if not os.path.isfile(lookin + '/mpc.h'):
                use_mpc = False
            break
    for adir in dirord:
        lookin = '%s/lib' % adir
        if os.path.isfile(lookin + '/lib' + mplib + '.a'):
            libdirs = [lookin]
            if local_dir and lookin.startswith(local_dir):
                rundirs = [lookin]
            # Verify that MPFR and MPC exist in the same directory
            if not os.path.isfile(lookin + '/libmpfr.a'):
                use_mpfr = False
            if not os.path.isfile(lookin + '/libmpc.a'):
                use_mpc = False
            break

# Validate include and library directories on Windows
if running_build and sys.version.find('MSC') != -1:
    rundirs = None
    my_extra_link_args = ["/MANIFEST"]
    if not local_dir:
        sys.stdout.write("Please specify parent directory of include and lib\n")
        sys.stdout.write("directories using -Ddir=<path>.");
        sys.exit()
    testpath = local_dir + '\\include'
    if os.path.isfile(testpath + '\\' + mplib + '.h'):
        incdirs = [testpath]
        if not os.path.isfile(testpath + '\\mpfr.h'):
            use_mpfr = False
        if not os.path.isfile(testpath + '\\mpc.h'):
            use_mpc = False
    testpath = local_dir + '\\lib'
    if os.path.isfile(testpath + '\\' + mplib + '.lib'):
        libdirs = [testpath]
        if not os.path.isfile(testpath + '\\mpfr.lib'):
            use_mpfr = False
        if not os.path.isfile(testpath + '\\mpc.lib'):
            use_mpc = False

# Use options to prevent use of MPFR and MPC even if found.
#   -DNOMPC  -> build without MPC library
#   -DNOMPFR -> build without MPFR library
for token in sys.argv:
    if token.upper() == '-DNOMPC':
        use_mpc = False
    if token.upper() == '-DNOMPFR':
        use_mpfr = False

if not use_mpfr:
    use_mpc = False

# Configure the defines...
defines = []
if mplib == 'mpir':
    defines.append( ('MPIR', 1) )
if use_mpfr:
    defines.append( ('WITHMPFR', 1) )
else:
    defines.append( ('NOMPFR', 1) )
if use_mpc:
    defines.append( ('WITHMPC', 1) )
else:
    defines.append( ('NOMPC', 1) )

# Build list of the required libraries...
libs = [mplib]
if use_mpfr:
    libs.append('mpfr')
if use_mpc:
    libs.append('mpc')

# Error message if libraries can not be found...
if running_build and not libdirs:
    sys.stdout.write("GMPY2 can not find the requires libraries.\n")
    sys.exit()

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'
# prepare the extension for building
gmpy2_ext = Extension('gmpy2', sources=['src/gmpy2.c'],
    include_dirs=incdirs,
    library_dirs=libdirs,
    libraries=libs,
    runtime_library_dirs=rundirs,
    define_macros = defines,
    extra_link_args = my_extra_link_args)

setup (name = "gmpy2",
       version = "2.0.0b4",
       maintainer = "Case Van Horsen",
       maintainer_email = "casevh@gmail.com",
       url = "http://code.google.com/p/gmpy/",
       description = "GMP/MPIR, MPFR, and MPC interface to Python 2.6+ and 3.x",

       classifiers = [
         'Development Status :: 4 - Beta',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
         'Natural Language :: English',
         'Operating System :: MacOS :: MacOS X',
         'Operating System :: Microsoft :: Windows',
         'Operating System :: POSIX',
         'Programming Language :: C',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 3',
         'Topic :: Scientific/Engineering :: Mathematics',
         'Topic :: Software Development :: Libraries :: Python Modules',
       ],

       ext_modules = [ gmpy2_ext ]
)
