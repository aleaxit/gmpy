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

Changes in gmpy2 2.0.4
----------------------

* Fixed bit_scan0() for negative values.
* Added option to setup.py (--static) to support static linking.
* Manpage is now installed in section 3.

Changes in gmpy2 2.0.3
----------------------

* Fixed bugs in lucas2() and atanh() that caused incorrect results.

Changes in gmpy2 2.0.2
----------------------

* Rebuild the Windows binary installers due to a bug in MPIR.
* Correct test in is_extra_strong_lucas_prp(). Note: The incorrect test is not
  known to cause any errors.

Changes in gmpy2 2.0.1
----------------------

* Updated setup.py to work in more situations.
* Corrected exception handling in basic operations with mpfr type.
* Correct InvalidOperation exception not raised in certain circumstances.
* invert() now raises an exception if the modular inverse does not exist.
* Fixed internal exception in is_bpsw_prp() and is_strong_bpsw_prp().
* Updated is_extra_strong_lucas_prp() to latest version.

Changes in gmpy2 2.0.0
----------------------

* Fix segmentation fault in _mpmath_normalize (an undocumented helper function
  specifically for mpmath.)
* Improved setup.py See below for documentation on the changes.
* Fix issues when compiled without support for MPFR.
* Conversion of too large an mpz to float now raises OverflowError instead of
  returning *inf*.
* Renamed min2()/max2() to minnum()/maxnum()
* The build and install process (i.e. setup.py) has been completely rewritten.
  See the Installation section for more information.
* get_context() no longer accepts keyword arguments.

Known issues in gmpy2 2.0.0
-----------------------------

* The test suite is still incomplete.

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

Installation
============

Installing gmpy2 on Windows
---------------------------

Pre-compiled versions of gmpy2 are available at `PyPi
<https://pypi.python.org/pypi/gmpy2/>`_ . Please select the installer
that corresponds to the version of Python installed on your computer.
Note that either a 32 or 64-bit version of Python can be installed on a
64-bit version of Windows. If you get an error message stating that
Python could not be found in the registry, you have the wrong
version of the gmpy2 installer.

Installing gmpy2 on Unix/Linux
------------------------------

Requirements
^^^^^^^^^^^^

gmpy2 has only been tested with recent versions of GMP, MPFR and MPC.
Specifically, for integer and rational support, gmpy2 requires GMP 5.1.x or
later. To support multiple-precision floating point arithmetic, MPFR 3.1.x or
later is required. MPC 1.0.1 or later is required for complex arithmetic.

Short Instructions
^^^^^^^^^^^^^^^^^^

You will need to install the development libraries for Python, GMP, MPFR, and
MPC. Different Linux distributions may the development packages differently.
Typical names are libpython-dev, libgmp-dev, libmpfr-dev, and libmpc-dev.

If your system includes recent versions of GMP, MPFR and MPC, and you have
the development libraries installed, compiling should be as simple as:

::

    cd <gmpy2 source directory>
    python setup.py build
    sudo python setup.py install

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
your home directory. The examples will use /home/<username>/local. If you use
one of standard directories (say /opt/local), then you won't need to specify
--prefix=/home/case/local to setup.py but you will need to specify the prefix
when compiling GMP, MPFR, and MPC.

Please substitute your actual user name for <username>.

Create the desired destination directory for GMP, MPFR, and MPC.
::

    $ mkdir /home/<username>/local

Download and un-tar the GMP source code. Change to the GMP source directory and
compile GMP.
::

    $ cd /home/<username>/local/src/gmp-6.0.0
    $ ./configure --prefix=/home/<username>/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPFR source code. Change to the MPFR source directory
and compile MPFR.
::

    $ cd /home/<username>/local/src/mpfr-3.1.2
    $ ./configure --prefix=/home/<username>/local --with-gmp=/home/<username>/local
    $ make
    $ make check
    $ make install

Download and un-tar the MPC source code. Change to the MPC source directory
and compile MPC.
::

    $ cd /home/<username>/local/src/mpc-1.0.2
    $ ./configure --prefix=/home/<username>/local --with-gmp=/home/<username>/local --with-mpfr=/home/<username>/local
    $ make
    $ make check
    $ make install

Compile gmpy2 and specify the location of GMP, MPFR and MPC. The location of
the GMP, MPFR, and MPC libraries is embedded into the gmpy2 library so the new
versions of GMP, MPFR, and MPC do not need to be installed the system library
directories. The prefix directory is added to the beginning of the directories
that are checked so it will be found first.
::

    $ python setup.py install --prefix=/home/<username>/local

If you get a "permission denied" error message, you may need to use::

    $ python setup.py build --prefix=/home/<username>/local
    $ sudo python setup.py install --prefix=/home/<username>/local

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

**--prefix=<...>**
    Specify the directory prefix where GMP/MPIR, MPFR, and MPC are located. For
    example, **--prefix=/opt/local** instructs setup.py to search /opt/local/include
    for header files and /opt/local/lib for libraries.

**--nompfr**
    Disables support for MPFR and MPC. This option is intended for testing purposes
    and is not offically supported. MPFR will be required for future versions of
    gmpy2.

**--nompc**
    Disables support for MPC. This option is intended for testing purposes and is not
    officially supported. MPC will be required for future versions of gmpy2.

**--static**
    Builds a statically linked library. This option will likely require the use of
    --prefix=<...> to specify the directory where the static libraries are located.
    To successfully link with gmpy2, the GMP, MPFR, and MPC libraries must be compiled
    with the --with-pic option.

