r'''
>>> gmpy.version()
'1.15'
>>>
'''

import sys
import doctest
import gmpy

if sys.argv[-1] == 'debug':
    gmpy.set_debug(1)

import gmpy_test_cvr
import gmpy_test_rnd
import gmpy_test_mpf
import gmpy_test_mpq
import gmpy_test_mpz
import gmpy_test_dec


test_modules = (gmpy_test_cvr, gmpy_test_rnd, gmpy_test_mpf,
    gmpy_test_mpq, gmpy_test_mpz, gmpy_test_dec)

_g = gmpy
print("Unit tests for gmpy 1.15")
print("    on Python %s" % sys.version)
if _g.gmp_version():
    print("Testing gmpy %s (GMP %s), default caching (%s, %s)" % (
        (_g.version(), _g.gmp_version(), _g.get_cache()[0],
        _g.get_cache()[1])))
else:
    print("Testing gmpy %s (MPIR %s), default caching (%s, %s)" % (
        (_g.version(), _g.mpir_version(), _g.get_cache()[0],
        _g.get_cache()[1])))

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
    print("Python's Fraction and gmpy's mpq.")
