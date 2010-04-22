import sys, os
from distutils.core import setup, Extension

# Fail gracefully for old versions of Python.
if sys.version < '2.6.0':
    sys.stdout.write("GMPY2 requires Python 2.6 or later.\n")
    sys.stdout.write("Please use GMPY 1.x for earlier versions of Python.\n")
    sys.exit()

# Check if MPIR or GMP should be used.
mplib='gmp'
for token in sys.argv:
    if token.startswith('-D') and 'MPIR' in token:
        mplib='mpir'
        break

# determine include and library dirs
incdirs = libdirs = ()
if sys.version.find('MSC') == -1:
    # Unix-like build (including MacOSX)
    incdirs = ['./src']
    dirord = ['/opt/local', '/usr/local']
    for adir in dirord:
        lookin = '%s/include' % adir
        if os.path.isfile(lookin + '/' + mplib + '.h'):
            incdirs.append(lookin)
            break
    for adir in dirord:
        lookin = '%s/lib' % adir
        if os.path.isfile(lookin + '/lib' + mplib + '.a'):
            libdirs = [lookin]
            break

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'
# prepare the extension for building
gmpy2_ext = Extension('gmpy2', sources=['src/gmpy2.c'],
    include_dirs=incdirs,
    library_dirs=libdirs,
    libraries=[mplib])

setup (name = "gmpy2",
       version = "2.0.0a0",
       maintainer = "Alex Martelli",
       maintainer_email = "aleaxit@gmail.com",
       url = "http://code.google.com/p/gmpy/",
       description = "GMP/MPIR interface to Python 2.5+ and 3.x",

       classifiers = [
         'Development Status :: 5 - Production/Stable',
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
