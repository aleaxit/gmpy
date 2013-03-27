import sys
import doctest
import gmpy2

if sys.argv[-1] == 'debug':
    gmpy.set_debug(1)

import gmpy_test_cvr
#~ import gmpy_test_rnd
import gmpy_test_mpf
import gmpy_test_mpq
import gmpy_test_mpz
import gmpy_test_dec
import gmpy_test_xmpz


test_modules = (gmpy_test_cvr, gmpy_test_mpf,
    gmpy_test_mpq, gmpy_test_mpz, gmpy_test_dec, gmpy_test_xmpz)

_g = gmpy2
print("Unit tests for gmpy2")
print("    on Python %s" % sys.version)
print("Testing gmpy2 {0}".format(_g.version()))
print("  Mutliple-precision library:   {0}".format(_g.mp_version()))
print("  Floating-point library:       {0}".format(_g.mpfr_version()))
print("  Complex library:              {0}".format(_g.mpc_version()))
print("  Caching Values: (Number)      {0}".format(_g.get_cache()[0]))
print("  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1]))

pf, pt = 0, 0
for x in test_modules:
    testit = x._test()
    failures, tests = testit
    if tests == 0: continue
    print("%s %3d tests, %d failures" % (x.__name__, tests-pt, failures-pf))
    pf, pt = failures, tests

doctest.master.summarize(1)

if sys.version_info < (3,1,1):
    print("There is a known bug with Fraction == mpq for versions of Python")
    print("less than 3.1.1. Please upgrade if you rely on comparisons between")
    print("Python's Fraction and gmpy2's mpq.")
