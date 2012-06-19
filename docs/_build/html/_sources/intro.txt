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

Installing gmpy2 on Windows
---------------------------

Pre-compiled versions of gmpy2 are available at `Downloads
<http://code.google.com/p/gmpy/downloads/list>`_ . Please
select the installer that corresponds to the version of Python installed on
your computer. Note that either a 32 or 64-bit version of Python can be
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
later is required. MPC 0.9 or later is required for complex arithmetic.

The MPC and MPFR libraries are optional. If the MPC library is not available,
gmpy2 will still support integer, rational, and real floating-point arithmetic.
If the MPFR library is not available, gmpy2 will only support integer and
rational arithmetic. The mpf type included with GMP is no longer supported.

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
with existing libraries on your system, the following instructions install GMP,
MPFR and MPC in a separate directory. The examples use /opt/local but you can
use another directory if you choose.

Create the desired destination directory for GMP, MPFR, and MPC.
::

$ mkdir /opt/local

Download and un-tar the GMP source code. Change to GMP source directory and
compile GMP.
::

    $ cd /opt/local/src/gmp-5.0.2
    $ ./configure --prefix=/opt/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPFR source code. Change to MPFR source directory
and compile MPFR.
::

    $ cd /opt/local/mpfr-3.1.0
    $ ./configure --prefix=/opt/local --with-gmp=/opt/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPC source code. Change to MPC source directory
and compile MPC.
::

    $ cd /opt/local/mpc-0.9
    $ ./configure --prefix=/opt/local --with-gmp=/opt/local --with-mpfr=/opt/local
    $ make
    $ make check
    $ make install

Compile gmpy2 and specify the location of GMP, MPFR and MPC.
::

    $ python setup.py build_ext -Ddir=/opt/local install

If you get a "permission denied" error message, you may need to use::

    $ sudo python setup.py build_ext -Ddir=/home/opt/local install

Miscellaneous gmpy2 Functions
-----------------------------

**from_binary(...)**
    from_binary(bytes) returns a gmpy2 object from a byte sequence created by
    to_binary().

**get_cache(...)**
    get_cache() returns the current cache size (number of objects) and the
    maximum size per object (number of limbs).

    gmpy2 maintains an internal list of freed *mpz*, *xmpz*, *mpq*, and *mpfr*
    objects for reuse. The cache significantly improves performance but does
    increase the memory footprint.

**license(...)**
    license() returns the gmpy2 license information.

**mp_limbsize(...)**
    mp_limbsize() returns the number of bits per limb used by the GMP or MPIR
    libarary.

**mp_version(...)**
    mp_version() returns the version of the GMP or MPIR library.

**mpc_version(...)**
    mpc_version() returns the version of the MPC library.

**mpfr_version(...)**
    mpfr_version() returns the version of the MPFR library.

**random_state(...)**
    random_state([seed]) returns a new object containing state information for
    the random number generator. An optional integer argument can be specified
    as the seed value. Only the Mersenne Twister random number generator is
    supported.

**set_cache(...)**
    set_cache(number, size) updates the maximum number of freed objects of each
    type that are cached and the maximum size (in limbs) of each object. The
    maximum number of objects of each type that can be cached is 1000. The
    maximum size of an object is 16384. The maximum size of an object is
    approximately 64K on 32-bit systems and 128K on 64-bit systems.

    .. note::
        The caching options are global to gmpy2. Changes are not thread-safe. A
        change in one thread will impact all threads.

**to_binary(...)**
    to_binary(x) returns a byte sequence from a gmpy2 object. All object types
    are supported.

**version(...)**
    version() returns the version of gmpy2.
