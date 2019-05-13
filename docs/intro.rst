Introduction to gmpy2
=====================

gmpy2 is a C-coded Python extension module that supports multiple-precision
arithmetic. gmpy2 is the successor to the original gmpy module. The gmpy module
only supported the GMP multiple-precision library. gmpy2 adds support for the
MPFR (correctly rounded real floating-point arithmetic) and MPC (correctly
rounded complex floating-point arithmetic) libraries. gmpy2 also updates the
API and naming conventions to be more consistent and support the additional
functionality.

The following libraries are supported:

* GMP for integer and rational arithmetic

  Home page: http://gmplib.org
* MPIR is based on the GMP library but adds support for Microsoft's Visual
  Studio compiler. It is used to create the Windows binaries.

  Home page: http://www.mpir.org
* MPFR for correctly rounded real floating-point arithmetic

  Home page: http://www.mpfr.org
* MPC for correctly rounded complex floating-point arithmetic

  Home page: http://mpc.multiprecision.org
* Generalized Lucas sequences and primality tests are based on the following
  code:

  mpz_lucas: http://sourceforge.net/projects/mpzlucas/

  mpz_prp: http://sourceforge.net/projects/mpzprp/

gmpy2 Versions
--------------

This manual documents the two major versions of gmpy2. Sections that are
specific to a particular version will be identified as such.

The 2.0 version is the stable release that only receives bug fixes and very
minor updates. Version 2.1 is currently under active development and includes
several new capabilities. Most gmpy2 2.0 code should run unchanged with
gmpy2 2.1.

The most significant change in gmpy2 2.1 is support for thread-safe contexts.
This change required extensive refactoring of almost all internal functions.

Please see the History chapter for a detail list of the changes.

Installation
============

Installing gmpy2 on Windows
---------------------------

Pre-compiled versions of gmpy2 2.0.8 are available at
`https://pypi.org/project/gmpy2/`.

A pre-compiled version of gmpy2 2.1.0a1 is available at
`https://pypi.org/project/gmpy2/2.1.0a1/`. Updated Windows versions should be
available again beginning with version 2.1.0b1.

Installing gmpy2 on Unix/Linux
------------------------------

Requirements
^^^^^^^^^^^^

gmpy2 has only been tested with the most recent versions of GMP, MPFR and MPC.
Specifically, for integer and rational support, gmpy2 requires GMP 5.0.x or
later. To support multiple-precision floating point arithmetic, MPFR 3.1.x or
later is required. MPC 1.0.1 or later is required for complex arithmetic.

Short Instructions
^^^^^^^^^^^^^^^^^^

gmpy2 requires the development files for GMP, MPFR, and MPC. The actual package
that provides these files varies between Linux distributions. Installing
"libmpc-dev" (or its equivalent) is usually sufficient.

If your system has the development libraries installed, compiling should be as
simple as:

::

    cd <gmpy2 source directory>
    python setup.py build_ext --force install --force

If this fails, read on.

Detailed Instructions
^^^^^^^^^^^^^^^^^^^^^

Note: You really shouldn't need to do this. Unless you need the capabilities
provided a newer GMP/MPFR/MPC, you should use the versions provided by your
distribution.

If your Linux distribution does not support recent versions of GMP, MPFR and
MPC, you will need to compile your own versions. To avoid any possible conflict
with existing libraries on your system, it is recommended to use a directory
not normally used by your distribution.

Create the desired destination directory for GMP, MPFR, and MPC.
::

    $ mkdir /home/<<your username>>/local

Download and un-tar the GMP source code. Change to the GMP source directory and
compile GMP.
::

    $ cd /home/<<your username>>/local/src/gmp-6.1.2
    $ ./configure --prefix=/home/<<your username>>/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPFR source code. Change to the MPFR source directory
and compile MPFR.
::

    $ cd /home/<<your username>>/local/src/mpfr-4.0.1
    $ ./configure --prefix=/home/<<your username>>/local --with-gmp=/home/<<your username>>/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPC source code. Change to the MPC source directory
and compile MPC.
::

    $ cd /home/<<your username>>/local/src/mpc-1.1.0
    $ ./configure --prefix=/home/<<your username>>/local --with-gmp=/home/<<your username>>/local --with-mpfr=/home/<<your username>>/local
    $ make
    $ make check
    $ make install

Compile gmpy2 and specify the location of GMP, MPFR and MPC. The location of
the GMP, MPFR, and MPC libraries is embedded into the gmpy2 library so the new
versions of GMP, MPFR, and MPC do not need to be installed the system library
directories. The prefix directory is added to the beginning of the directories
that are checked so it will be found first.
::

    $ python setup.py install --prefix=/home/case/local

If you get a "permission denied" error message, you may need to use::

    $ python setup.py build --prefix=/home/case/local
    $ sudo python setup.py install --prefix=/home/case/local

Options for setup.py
^^^^^^^^^^^^^^^^^^^^

**--force**
    Ignore the timestamps on all files and recompile. Normally, the results of a
    previous compile are cached. To force gmpy2 to recognize external changes
    (updated version of GMP, etc.), you will need to use this option.

**--mpir**
    Force the use of MPIR instead of GMP. GMP is the default library on non-Windows
    operating systems.

**--gmp**
    Force the use of GMP instead of MPIR. MPIR is the default library on Windows
    operating systems.

**--shared=<...>**
    Add the specified directory prefix to the beginning of the list of
    directories that are searched for GMP, MPFR, and MPC shared libraries.

**--static=<...>**
    Create a statically linked library using libraries from the specified path,
    or from the operating system's default library location if no path is specified

