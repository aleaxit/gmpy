from __future__ import print_function, division

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

# If mpc version < 1.1.0 gmpy2.root_of_unity is defined and gmpy2.cmp_abs
# doesn't manage complex parameters.
# We create a doctest flag to skip a doctest when mpc version is < 1.1.0
SKIP_MPC_LESS_THAN_110 = doctest.register_optionflag("SKIP_MPC_LESS_THAN_110")
mpc_version_110 = 'root_of_unity' in dir(gmpy2) # True if mpc version >= 1.1.0

SKIP_IN_DEBUG_MODE = doctest.register_optionflag("SKIP_IN_DEBUG_MODE")

class Gmpy2DocTestParser(DocTestParser):
    def parse(self, *args, **kwargs):
        examples = DocTestParser.parse(self, *args, **kwargs)
        for example in examples:
            if not isinstance(example, Example):
                continue
            if not mpc_version_110 and SKIP_MPC_LESS_THAN_110 in example.options:
                example.options[SKIP] = True
            if debug and SKIP_IN_DEBUG_MODE in example.options:
                example.options[SKIP] = True

        return examples

parser = Gmpy2DocTestParser()

print()
print("Unit tests for gmpy2 {0} with Python {1}".format(gmpy2.version(), sys.version.split()[0]))
print("  Mutliple-precision library:     {0}".format(gmpy2.mp_version()))
print("  Floating-point library:         {0}".format(gmpy2.mpfr_version()))
print("  Complex library:                {0}".format(gmpy2.mpc_version()))
print("  Caching Values: (Cache size)    {0}".format(gmpy2.get_cache()[0]))
print("  Caching Values: (Size in limbs) {0}".format(gmpy2.get_cache()[1]))
print()

mpz_doctests = ["test_mpz_create.txt", "test_mpz.txt", "test_mpz_io.txt",
                "test_mpz_pack_unpack.txt", "test_misc.txt"]

mpq_doctests = ["test_mpq.txt"]

mpfr_doctests = ["test_mpfr_create.txt", "test_mpfr.txt",
                 "test_mpfr_trig.txt", "test_mpfr_min_max.txt",
                 "test_context.txt", "test_mpfr_subnormalize.txt"]

# Some tests may differ between MPFR3 and MPFR4.
mpfr_major_version = gmpy2.mpfr_version().split()[1].split('.')[0]
mpfr_version_tests = [os.path.basename(i)
                      for i in glob.glob(os.path.join(os.path.dirname(__file__),
                                         "test_mpfr" + mpfr_major_version + "*.txt"))]

mpc_doctests = ["test_mpc_create.txt", "test_mpc.txt", "test_mpc_trig.txt"]

gmpy2_tests = [os.path.basename(i)
               for i in glob.glob(os.path.join(os.path.dirname(__file__),
                                  "test_gmpy2*.txt"))]

# The following tests will only pass on Python 3.2+.
py32_doctests = ["test_py32_hash.txt"]

failed = 0
attempted = 0

all_doctests = gmpy2_tests + mpz_doctests + mpq_doctests

all_doctests += mpfr_doctests + mpfr_version_tests

all_doctests += mpc_doctests

if sys.version_info > (3,1):
    all_doctests += py32_doctests

for test in sorted(all_doctests):
    if test.endswith("py2.txt") and sys.version_info[0] >= 3:
        continue
    if test.endswith("py3.txt") and sys.version_info[0] < 3:
        continue
    for r in range(repeat):
        result = doctest.testfile(test, globs=globals(),
                                  optionflags=doctest.IGNORE_EXCEPTION_DETAIL |
                                              doctest.NORMALIZE_WHITESPACE |
                                              doctest.REPORT_NDIFF,
                                  parser=parser)
        print("Results for:  {0:25}".format(test.split(".")[0]), end="")
        print(" Attempted: {1:4d}   Failed: {0:4d}".format(*result), end="")
        if debug:
            print(" RefCount: {0:6d}".format(sys.gettotalrefcount()))
        else:
            print()
        failed += result[0]
        attempted += result[1]
    if repeat > 1:
        print()


print()
print("                              Summary - Attempted: {0:4d}   Failed: {1:4d}".format(attempted, failed))
print()
print("Running external test programs.")

print("Running {0:30}  ".format("test_pack.py"), end="")
import test_pack
if test_pack.test():
    print("successful")
    attempted += 1
else:
    print("failed")
    failed += 1

print("Running {0:30}  ".format("test_mpz_args.py"), end="")
import test_mpz_args
if test_mpz_args.test():
    print("successful")
    attempted += 1
else:
    print("failed")
    failed += 1

if failed:
    sys.exit(1)
else:
    sys.exit(0)
