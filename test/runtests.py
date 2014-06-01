from __future__ import print_function, division

import sys
import doctest
import gmpy2

# *****************************************************************************
# Test strategy
# -------------
# Tests are divided into two different categories:
#
#   1) The 'txt' files contain doctest style tests. These tests should cover
#      basic functionality for all functions/types.
#   2) The 'py' files contain Python code that perform extensive tests, but
#      may not test every function.
#
# If run by a debug build of Python, the test suite can be repeated multiple
# times to search for memory leaks.
#
# NOTE: IF THE LAST TEST IN A BLOCK OF TESTS GENERATES AN EXCEPTION, THE
#       REFERENCE COUNTING IN A DEBUG BUILD GETS CONFUSED. ALWAYS ENSURE THAT
#       AT LEAST ONE VALID TEST IS PERFORMED AFTER AN EXCEPTION IS RAISED!
#
# *****************************************************************************

# Check if this is a debug build of Python.
try:
    sys.gettotalrefcount()
    debug = True
except AttributeError:
    debug = False

# Change repeat to the number of times to repeat each test. Combined with a
# debug build, this can help identify memory leaks.
if debug:
    try:
        repeat = abs(int(sys.argv[1]))
    except:
        repeat = 1
else:
    repeat = 1

print()
print("Unit tests for gmpy2 {0} with Python {1}".format(gmpy2.version(), sys.version.split()[0]))
print("  Mutliple-precision library:     {0}".format(gmpy2.mp_version()))
print("  Floating-point library:         {0}".format(gmpy2.mpfr_version()))
print("  Complex library:                {0}".format(gmpy2.mpc_version()))
print("  Caching Values: (Cache size)    {0}".format(gmpy2.get_cache()[0]))
print("  Caching Values: (Size in limbs) {0}".format(gmpy2.get_cache()[1]))
print()

if sys.version.startswith('3.1'):
    print("Due to differences in formatting of exceptions and Python 3.x, there")
    print("will be test failures for exception handling when the tests are run")
    print("with Python 3.1. The doctest module in Python 3.2 and later does not")
    print("have this issue.")
    print()
    input("Press ENTER to continue.. ")
    print()

# The following tests should pass on all builds.
mpz_doctests = ["test_mpz.txt", "test_mpz_io.txt", "test_mpz_pack_unpack.txt",
                "test_mpz_to_from_binary.txt"]
mpq_doctests = ["test_mpq.txt", "test_mpq_to_from_binary.txt"]

# The following tests require MPFR support.
mpfr_doctests = ["test_mpfr.txt", "test_mpfr_trig.txt", "test_mpfr_min_max.txt",
                 "test_mpfr_to_from_binary.txt", "test_context.txt"]

# The following tests require MPC support.
mpc_doctests = ["test_mpc.txt", "test_mpc_to_from_binary.txt"]

# The following tests will only pass on Python 3.2+.
py32_doctests = ["test_py32_hash.txt"]

failed = 0
attempted = 0

all_doctests = ["test_misc.txt"] + mpz_doctests + mpq_doctests
if gmpy2.mpfr_version():
    all_doctests += mpfr_doctests
if gmpy2.mpc_version():
    all_doctests += mpc_doctests
if sys.version >= "3.2":
    all_doctests += py32_doctests

for test in sorted(all_doctests):
    for r in range(repeat):
        result = doctest.testfile(test, globs=globals(), optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
        print("Results for:  {0:24}".format(test.split(".")[0]), end="")
        print(" Attempted: {1:4d}   Failed: {0:4d}".format(*result), end="")
        if debug:
            print(" RefCount: {0:6d}".format(sys.gettotalrefcount()))
        else:
            print()
        failed += result[0]
        attempted += result[1]


print()
print("                             Summary - Attempted: {0:4d}   Failed: {1:4d}".format(attempted, failed))

if failed:
    sys.exit(1)
else:
    sys.exit(0)

