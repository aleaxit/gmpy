from __future__ import print_function, division

import sys
import doctest
import gmpy2

print()
print("Unit tests for gmpy2 {0} with Python {1}".format(gmpy2.version(), sys.version.split()[0]))
print("  Mutliple-precision library:   {0}".format(gmpy2.mp_version()))
print("  Floating-point library:       {0}".format(gmpy2.mpfr_version()))
print("  Complex library:              {0}".format(gmpy2.mpc_version()))
print("  Caching Values: (Number)      {0}".format(gmpy2.get_cache()[0]))
print("  Caching Values: (Size, limbs) {0}".format(gmpy2.get_cache()[1]))
print()

# The following tests should pass on all supported versions of Python.
all_doctests = ["test_misc.txt", "test_dec.txt", "test_mpz.txt"]

# The following tests will only pass on Python 3.2+.
py32_doctests = ["test_hash.txt"]

failed = 0
attempted = 0

for test in all_doctests:
    print("Running test: {0}".format(test))
    result = doctest.testfile(test)
    print("  Results - Attempted: {1:4d} Failed: {0:4d}".format(*result))
    failed += result[0]
    attempted += result[1]

if sys.version_info >= (3,2):
    for test in py32_doctests:
        print("Running test: {0}".format(test))
        result = doctest.testfile(test)
        print("  Results - Attempted: {1:4d} Failed: {0:4d}".format(*result))
        failed += result[0]
        attempted += result[1]

print()
print("Summary of all tests")
print("  Summary - Attempted: {0:4d} Failed: {1:4d}".format(attempted, failed))

