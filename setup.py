import platform
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


_comp_args = ["-DSHARED=1"]
link_args = []
ON_WINDOWS = platform.system() == 'Windows'

sources = ['src/gmpy2.c']
_libs = ['gmp', 'mpfr', 'mpc'] if not ON_WINDOWS else [
    'libmpir', 'libmpfr', 'libmpc']


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
         "Windows"),
        ('static', None, "Enable static linking compile time options."
         "Static options are enabled by default on Windows platforms only."),
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
            _comp_args.append('-DFAST=1')
        if self.gcov:
            if ON_WINDOWS:
                raise ValueError("Cannot enable GCC code coverage on Windows")
            _comp_args.append('-DGCOV=1')
            _comp_args.append('-O0')
            _comp_args.append('--coverage')
            link_args.append('--coverage')
            _libs.append('gcov')
        if self.vector:
            _comp_args.append('-DVECTOR=1')
        if self.mpir or ON_WINDOWS:
            _comp_args.append('-DMPIR=1')
        if self.mpir and not ON_WINDOWS:
            _libs.remove('gmp')
            _libs.append('mpir')
        if self.static or ON_WINDOWS:
            _comp_args.remove('-DSHARED=1')
            _comp_args.append('-DSTATIC=1')
        if self.gdb:
            _comp_args.append('-ggdb')

    def build_extensions(self):
        compiler = self.compiler.compiler_type
        if compiler == 'mingw32':
            _comp_args.append('-DMSYS2=1')
        build_ext.build_extensions(self)


extensions = [
    Extension('gmpy2.gmpy2',
              sources=sources,
              include_dirs=["/usr/local/include", './src'],
              libraries=_libs,
              extra_compile_args=_comp_args,
              extra_link_args=link_args,
              )
]

cmdclass = {'build_ext': Gmpy2Build}

setup(
    name="gmpy2",
    version="2.1.0a1",
    author="Case Van Horsen",
    author_email="casevh@gmail.com",
    cmdclass=cmdclass,
    license="LGPL-3.0+",
    url="https://github.com/aleaxit/gmpy",
    description="gmpy2 interface to GMP/MPIR, MPFR, "
    "and MPC for Python 2.6+ and 3.4+",
    zip_safe=False,
    include_package_data=True,
    package_data={'gmpy2': [
        'gmpy2.pxd',
        'gmpy2.h',
    ]},
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
