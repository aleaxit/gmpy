r'''
>>> gmpy.version()
'1.15'
>>>
'''

import sys
import doctest
import gmpy

def writeln(s):
    sys.stdout.write(s+'\n')

if sys.version_info[0] == 3:
    writeln("Please use 'test3/gmpy_test.py' to test with Python 3.x.")
    sys.exit(0)

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
writeln("Unit tests for gmpy 1.15")
writeln("    on Python %s" % sys.version)
if _g.gmp_version():
    writeln("Testing gmpy %s (GMP %s), default caching (%s, %s)" % (
            (_g.version(), _g.gmp_version(), _g.get_cache()[0],
            _g.get_cache()[1])))
else:
    writeln("Testing gmpy %s (MPIR %s), default caching (%s, %s)" % (
            (_g.version(), _g.mpir_version(), _g.get_cache()[0],
            _g.get_cache()[1])))


pf, pt = 0, 0
for x in test_modules:
    testit = x._test()
    failures, tests = testit
    if tests == 0: continue
    writeln("%s %3d tests, %d failures" % (x.__name__, tests-pt, failures-pf))
    pf, pt = failures, tests

doctest.master.summarize(1)
