from __future__ import print_function, division

import sys
import doctest
import gmpy2

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

# The following tests should pass on all builds.
mpz_doctests = ["test_misc.txt"]

# The following tests require MPFR support.
mpfr_doctests = ["test_mpfr_trig.txt"]

# The following tests require MPC support.
mpc_doctests = []

# The following tests will only pass on Python 3.2+.
py32_doctests = ["test_py32_hash.txt"]

failed = 0
attempted = 0

all_doctests = mpz_doctests + mpfr_doctests + mpc_doctests
for test in all_doctests:
    print("Running test: {0:16}".format(test.split(".")[0]))
    for r in range(repeat):
        result = doctest.testfile(test, globs=globals())
        print("Results for:  {0:16}".format(test.split(".")[0]), end="")
        print(" Attempted: {1:4d} Failed: {0:4d}".format(*result), end="")
        if debug:
            print(" RefCount: {0:6d}".format(sys.gettotalrefcount()))
        else:
            print()
        failed += result[0]
        attempted += result[1]

if sys.version_info >= (3,2):
    for test in py32_doctests:
        print("Running test: {0:16}".format(test.split(".")[0]))
        for r in range(repeat):
            result = doctest.testfile(test, globs=globals())
            print("Results for:  {0:16}".format(test.split(".")[0]), end="")
            print(" Attempted: {1:4d} Failed: {0:4d}".format(*result), end="")
            if debug:
                print(" RefCount: {0:6d}".format(sys.gettotalrefcount()))
            else:
                print()
            failed += result[0]
            attempted += result[1]

print()
print("                     Summary - Attempted: {0:4d} Failed: {1:4d}".format(attempted, failed))

