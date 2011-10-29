from __future__ import print_function

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

doctests = ["test_misc.txt"]
failed = 0
attempted = 0

for test in doctests:
    print("Running test: {0}".format(test))
    result = doctest.testfile(test)
    print("  Results - Attempted: {1:4d} Failed: {0:4d}".format(*result))
    failed += result[0]
    attempted += result[1]

print()
print("Summary of all tests")
print("  Summary - Attempted: {0:4d} Failed: {1:4d}".format(attempted, failed))

