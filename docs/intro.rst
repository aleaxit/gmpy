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

Changes in gmpy2 2.0.0b4
------------------------

* Added __ceil__, __floor__, __trunc__, and __round__ methods to mpz and mpq
  types.
* Added __complex__ to mpc type.
* round(mpfr) now correctly returns an mpz type.
* If no arguments are given to mpz, mpq, mpfr, mpc, and xmpz, return 0 of the
  appropriate type.
* Fix broken comparison between mpz and mpq when mpz is on the left.
* Added __sizeof__ to all types. *Note: sys.getsizeof() calls __sizeof__ to get
  the memory size of a gmpy2 object. The returned value reflects the size of the
  allocated memory which may be larger than the actual minimum memory required
  by the object.*

Known issues in gmpy2 2.0.0b4
-----------------------------

* The new test suite (test/runtest.py) is incomplete and some tests fail on
  Python 2.x due to formating issues.


Changes in gmpy2 2.0.0b3
------------------------

* mp_version(), mpc_version(), and mpfr_version() now return normal strings on
  Python 2.x instead of Unicode strings.
* Faster conversion of the standard library Fraction type to mpq.
* Improved conversion of the Decimal type to mpfr.
* Consistently return OverflowError when converting "inf".
* Fix mpz.__format__() when the format code includes "#".
* Add is_infinite() and deprecate is_inf().
* Add is_finite() and deprecate is_number().
* Fixed the various is_XXX() tests when used with mpc.
* Added caching for mpc objects.
* Faster code path for basic operation is both operands are mpfr or mpc.
* Fix mpfr + float segmentation fault.

Changes in gmpy2 2.0.0b2
------------------------

* Allow xmpz slice assignment to increase length of xmpz instance by specifying
  a value for stop.
* Fixed reference counting bug in several is_xxx_prp() tests.
* Added iter_bits(), iter_clear(), iter_set() methods to xmpz.
* Added powmod() for easy access to three argument pow().
* Removed addmul() and submul() which were added in 2.0.0b1 since they are
  slower than just using Python code.
* Bug fix in gcd_ext when both arguments are not mpz.
* Added ieee() to create contexts for 32, 64, or 128 bit floats.
* Bug fix in context() not setting emax/emin correctly if they had been changed
  earlier.
* Contexts can be directly used in with statement without requiring
  set_context()/local_context() sequence.
* local_context() now accepts an optional context.

Changes in gmpy2 2.0.0b1 and earlier
------------------------------------

* Renamed functions that manipulate individual bits to bit_XXX() to align with
  bit_length().
* Added caching for mpq.
* Added rootrem(), fib2(), lucas(), lucas2().
* Support changed hash function in Python 3.2.
* Added is_even(), is_odd().
* Add caching of the calculated hash value.
* Add xmpz (mutable mpz) type.
* Fix mpq formatting issue.
* Add read/write bit access using slices to xmpz.
* Add read-only bit access using slices to mpz.
* Add pack()/unpack() methods to split/join an integer into n-bit chunks.
* Add support for MPFR (casevh)
* Removed fcoform float conversion modifier.
* Add support for MPC.
* Added context manager.
* Allow building with just GMP/MPIR if MPFR not available.
* Allow building with GMP/MPIR and MPFR if MPC not available.
* Removed most instance methods in favor of gmpy2.function. The general guideline
  is that *properties* of an instance can be done via instance methods but
  *functions* that return a new result are done using gmpy2.function.
* Added __ceil__, __floor__, and __trunc__ methods since they are called by
  math.ceil(), math.floor(), and math.trunc().
* Removed gmpy2.pow() to avoid conflicts.
* Removed gmpy2._copy and added xmpz.copy.
* Added support for __format__.
* Added as_integer_ratio, as_mantissa_exp, as_simple_fraction.
* Updated rich_compare.
* Require MPFR 3.1.0+ to get divby0 support.
* Added fsum(), degrees(), radians().
* Updated random number generation support.
* Changed license to LGPL 3+.
* Added lucasu, lucasu_mod, lucasv, and lucasv_mod.
  *Based on code contributed by David Cleaver.*
* Added probable-prime tests.
  *Based on code contributed by David Cleaver.*
* Added to_binary()/from_binary.
* Renamed numdigits() to num_digits().
* Added keyword precision to constants.
* Added addmul() and submul().
* Added __round__(), round2(), round_away() for mpfr.
* round() is no longer a module level function.
* Renamed module functions min()/max() to min2()/max2().
*    No longer conflicts with builtin min() and max()
* Removed set_debug() and related functionality.


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
later is required. MPC 1.0.1 or later is required for complex arithmetic.

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

    $ cd /opt/local/src/gmp-5.1.0
    $ ./configure --prefix=/opt/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPFR source code. Change to MPFR source directory
and compile MPFR.
::

    $ cd /opt/local/mpfr-3.1.1
    $ ./configure --prefix=/opt/local --with-gmp=/opt/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPC source code. Change to MPC source directory
and compile MPC.
::

    $ cd /opt/local/mpc-1.0.1
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

    gmpy2 maintains an internal list of freed *mpz*, *xmpz*, *mpq*, *mpfr*, and
    *mpc* objects for reuse. The cache significantly improves performance but
    also increases the memory footprint.

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
