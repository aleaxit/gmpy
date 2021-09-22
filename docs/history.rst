Changes for gmpy2 releases
==========================

Changes in gmpy2 2.1.0rc2
-------------------------

* Documentation updates.
* Improvements to build environment.

Changes in gmpy2 2.1.0rc1
-------------------------

* Added support for embedded underscore characters in string
      literals.
* Allow GIL release for mpz/xmpz/mpq types only.

Changes in gmpy2 2.1.0b6
------------------------

* Improve argument type processing by saving type information to 
      decrease the number of type check calls. Especially helpful
      for mpfr and mpc types. (Not complete but common operations
      are done.)
* Resolve bug in mpfr to mpq conversion; issue #287.
* Added limited support for releasing the GIL; disabled by default;
      see context.allow_release_gil.
* Refactored handling of inplace operations for mpz and xmpz types;
      inplace operations on xmpz will only return an xmpz result.
* Refactored handling of conversion to C integer types. Some
      exception types changes to reflect Python types.
* gcd() and lcm() now support more than two arguments to align with
      the corresponding functions in the math module.
 
Changes in gmpy2 2.1.0b5
------------------------

* Avoid MPFR bug in mfr_fac_ui (gmpy2.factorial) on platforms where
      long is 32-bits and argument is >= 44787929.
* Fixed testing bugs with Python 2.7.
* Fixed mpz(0) to C long or long long.
* Fixed incorrect results in f2q().
* Adjust test suite to reflect changes in output in MPFR 4.1.0.

Changes in gmpy2 2.1.0b4
------------------------

* Fix comparisons with mpq and custom rational objects.
* Fixes for some uncommon integer conversions scenarios.
 
Changes in gmpy2 2.1.0b3
------------------------

* Version bump only.
 
Changes in gmpy2 2.1.0b2
------------------------

* Many bug fixes.
 
Changes in gmpy2 2.1.0b1
------------------------

* Added cmp() and cmp_abs().
* Improved compatibility with _numbers_ protocol.
* Many bug fixes.

Changes in gmpy2 2.1.a05
------------------------

* Fix qdiv() not returning mpz() when it should.
* Added root_of_unity().

Changes in gmpy2 2.1.0a4
------------------------

* Fix issue 204; missing file for Cython.
* Additional support for MPFR 4
   - Add fmma() and fmms()

Changes in gmpy2 2.1.0a3
------------------------

* Updates to setup.py.
* Initial support for MPFR4
   - Add nrandom()
   - grandom() now calls nrandom twice; may return different values versus
     MPFR3
   - Add rootn(); same as root() except different sign when taking even root
     of -0.0

Changes in gmpy2 2.1.0a2
------------------------

* Revised build process.
* Removal of unused code/macros.
* Cleanup of Cython interface.

Changes in gmpy2 2.1.0a1
------------------------

* Thread-safe contexts are now supported. Properly integrating thread-safe
  contexts required an extensive rewrite of almost all internal functions.
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
* If a Python object has an __mpz__ method, it will be called bye *mpz()* to
  allow an unrecognized type to be converted to an mpz instance. The same is
  true for the other gmpy2 types.
* A new C-API and Cython interface has been added.

Changes in gmpy2 2.0.4
----------------------

* Fix bit_scan0() for negative values.
* Changes to setup.py to allow static linking.
* Fix performance regression with mpmath and Python 3.

Changes in gmpy2 2.0.3
----------------------

* Fix lucas2() and atanh(); they were returning incorrect values.

Changes in gmpy2 2.0.2
----------------------

* Rebuild Windows binary installers due to MPIR 2.6.0 bug in next_prime().
* Another fix for is_extra_strong_lucas_prp().

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
  Python 2.x due to formatting issues.


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

