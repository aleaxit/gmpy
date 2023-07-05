import sys
import os
import glob
import doctest
from doctest import DocTestParser, Example, SKIP
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
# *****************************************************************************

test_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
print()
print("Unit tests for gmpy2 {0} with Python {1}".format(gmpy2.version(), sys.version.split()[0]))
print("  Mutliple-precision library:     {0}".format(gmpy2.mp_version()))
print("  Floating-point library:         {0}".format(gmpy2.mpfr_version()))
print("  Complex library:                {0}".format(gmpy2.mpc_version()))
print()

mpz_doctests = ["test_mpz_create.txt", "test_mpz.txt", "test_mpz_io.txt"]

mpq_doctests = ["test_mpq.txt"]

mpfr_doctests = ["test_mpfr_create.txt", "test_mpfr.txt",
                 "test_mpfr_trig.txt",
                 "test_context.txt", "test_mpfr_subnormalize.txt"]

# Some tests may differ between MPFR3 and MPFR4.
mpfr_major_version = gmpy2.mpfr_version().split()[1].split('.')[0]
mpfr_version_tests = [os.path.basename(i)
                      for i in glob.glob(os.path.join(test_dir,
                                         "test_mpfr" + mpfr_major_version + "*.txt"))]

mpc_doctests = ["test_mpc.txt", "test_mpc_trig.txt"]

gmpy2_tests = [os.path.basename(i)
               for i in glob.glob(os.path.join(test_dir,
                                  "test_gmpy2*.txt"))]

failed = 0
attempted = 0

all_doctests = gmpy2_tests + mpz_doctests + mpq_doctests

all_doctests += mpfr_doctests + mpfr_version_tests

all_doctests += mpc_doctests

for test in sorted(all_doctests):
    result = doctest.testfile(test, globs=globals(),
                                optionflags=doctest.IGNORE_EXCEPTION_DETAIL |
                                            doctest.NORMALIZE_WHITESPACE |
                                            doctest.REPORT_NDIFF)
    print("Results for:  {0:25}".format(test.split(".")[0]), end="")
    print(" Attempted: {1:4d}   Failed: {0:4d}".format(*result), end="")
    print()
    failed += result[0]
    attempted += result[1]


print()
print("                              Summary - Attempted: {0:4d}   Failed: {1:4d}".format(attempted, failed))
print()
print("Running external test programs.")

print("Running {0:30}  ".format("pytest"), end="")
if os.system(sys.executable + " -m pytest " + " ".join(glob.glob(test_dir + "/test_*.py"))) == 0:
    print("successful")
    attempted += 1
else:
    print("failed")
    failed += 1

if failed:
    sys.exit(1)
else:
    sys.exit(0)
