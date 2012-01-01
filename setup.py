import sys, os
from distutils.core import setup, Extension

# monkey-patch distutils if it can't cope with the "classifiers" and
# "download_url" keywords
if sys.version < '2.2.3':
    from distutils.dist import DistributionMetadata
    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None

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
    dirord = ['/usr', '/opt/local', '/usr/local']
    for adir in dirord:
        lookin = '%s/include' % adir
        if os.path.isfile(lookin + '/' + mplib + '.h'):
            incdirs.append(lookin)
            dirord = [adir]
            break
    for adir in dirord:
        lookin = '%s/lib' % adir
        if os.path.isfile(lookin + '/lib' + mplib + '.a'):
            libdirs = [lookin]
            break

# decomment next line (w/gcc, only!) to support gcov
#   os.environ['CFLAGS'] = '-fprofile-arcs -ftest-coverage -O0'
# prepare the extension for building
gmpy_ext = Extension('gmpy', sources=['src/gmpy.c'],
    include_dirs=incdirs,
    library_dirs=libdirs,
    libraries=[mplib])

setup (name = "gmpy",
       version = "1.15",
       maintainer = "Case Van Horsen",
       maintainer_email = "casevh@gmail.com",
       url = "http://code.google.com/p/gmpy/",
       description = "MPIR/GMP interface to Python 2.4+ and 3.x",

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

       ext_modules = [ gmpy_ext ]
)
