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

There are two versions of gmpy2. The 2.0 version is the stable release that
only receives bug fixes and very minor updates. Version 2.1 is currently under
active development and includes several new capabilities. Most gmpy2 2.0 code
should run unchanged with gmpy2 2.1

Enhancements in gmpy2 2.1
-------------------------

The most significant changes in gmpy2 2.1 are:

Changes in gmpy2 2.1.0a1
------------------------

* Thread-safe contexts are now supported. Properly integrating thread-safe
  contexts required an extensive rewrite of almost all internal functions.
  Changing the active context in one thread will no longer change the behavior
  in other threads.
* MPFR and MPC are now required. It is no longer possible to build a version
  of gmpy2 that only supports the GMP library.
* The function inverse() now raises an exception if the inverse does not
  exist.
* Context methods have been added for MPFR/MPC related functions.
* A new context option (*rational_division*) has been added that changes the
  behavior of integer division involving *mpz* instances to return a rational
  result instead of a floating point result.
* gmpy2 types are now registered in the numeric tower.
* In previous versions of gmpy2, *gmpy2.mpz* was a factory function that
  returned an  *mpz* instance. *gmpy2.mpz* is now an actual type. The same
  is true for the other gmpy2 types.
* If a Python object has an __mpz__ method, it will be called by *mpz()* to
  allow an unrecognized type to be converted to an mpz instance. The same is
  true for the other gmpy2 types.
* Support for Cython via the addition of a C-API and a gmpy2.pxd file.

Installation
============

This section will be updated soon to reflect improved support of pip and the
wheel format.

Installing gmpy2 on Windows
---------------------------

Pre-compiled versions of gmpy2 are available at `https://pypi.python.org/pypi/gmpy2`.
Please select the installer that corresponds to the version of Python installed
on your computer. Note that either a 32 or 64-bit version of Python can be
installed on a 64-bit version of Windows. If you get an error message
stating that Python could not be found in the registry, you have the wrong
version of the gmpy2 installer.

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

If your system includes sufficiently recent versions of GMP, MPFR and MPC, and
you have the development libraries installed, compiling should be as simple as:

::

    cd <gmpy2 source directory>
    python setup.py install

If this fails, read on.

Detailed Instructions
^^^^^^^^^^^^^^^^^^^^^

If your Linux distribution does not support recent versions of GMP, MPFR and
MPC, you will need to compile your own versions. To avoid any possible conflict
with existing libraries on your system, it is recommended to use a directory
not normally used by your distribution. setup.py will automatically search the
following directories for the required libraries:

    #. /opt/local
    #. /opt
    #. /usr/local
    #. /usr
    #. /sw

If you can't use one of these directories, you can use a directory located in
your home directory. The examples will use /home/case/local. If you use one of
standard directories (say /opt/local), then you won't need to specify
--prefix=/home/case/local to setup.py but you will need to specify the prefix
when compiling GMP, MPFR, and MPC.

Create the desired destination directory for GMP, MPFR, and MPC.
::

    $ mkdir /home/case/local

Download and un-tar the GMP source code. Change to the GMP source directory and
compile GMP.
::

    $ cd /home/case/local/src/gmp-5.1.0
    $ ./configure --prefix=/home/case/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPFR source code. Change to the MPFR source directory
and compile MPFR.
::

    $ cd /home/case/local/src/mpfr-3.1.1
    $ ./configure --prefix=/home/case/local --with-gmp=/home/case/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPC source code. Change to the MPC source directory
and compile MPC.
::

    $ cd /home/case/local/src/mpc-1.0.1
    $ ./configure --prefix=/home/case/local --with-gmp=/home/case/local --with-mpfr=/home/case/local
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

