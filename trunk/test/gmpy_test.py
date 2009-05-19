r'''
>>> gmpy.version()
'1.05'
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
print "Unit tests for gmpy 1.05"
print "    on Python %s" % sys.version
print "Testing gmpy %s (GMP %s), default caching (%s, %s, %s..%s)" % (
    (_g.version(), _g.gmp_version(), _g.get_zcache(), _g.get_qcache(),
            ) + _g.get_zconst())

pf, pt = 0, 0
for x in test_modules:
    testit = x._test()
    failures, tests = testit
    if tests == 0: continue
    print x.__name__,
    print "%3d tests, %d failures" % (tests-pt, failures-pf)
    pf, pt = failures, tests

doctest.master.summarize(1)
