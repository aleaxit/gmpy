import sys
import doctest
import gmpy2

def writeln(s):
    sys.stdout.write(s+'\n')

if sys.version_info[0] == 3:
    writeln("Please use 'test3/gmpy_test.py' to test with Python 3.x.")
    sys.exit(0)

if sys.argv[-1] == 'debug':
    gmpy.set_debug(1)

import gmpy_test_cvr
#~ import gmpy_test_rnd
import gmpy_test_mpf
import gmpy_test_mpq
import gmpy_test_mpz
import gmpy_test_dec


test_modules = (gmpy_test_cvr, gmpy_test_mpf,
    gmpy_test_mpq, gmpy_test_mpz, gmpy_test_dec)

_g = gmpy2
writeln("Unit tests for gmpy2")
writeln("    on Python %s" % sys.version)
writeln("Testing gmpy2 {0}".format(_g.version()))
writeln("  Mutliple-precision library:   {0}".format(_g.mp_version()))
writeln("  Floating-point library:       {0}".format(_g.mpfr_version()))
writeln("  Complex library:              {0}".format(_g.mpc_version()))
writeln("  Caching Values: (Number)      {0}".format(_g.get_cache()[0]))
writeln("  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1]))

pf, pt = 0, 0
for x in test_modules:
    testit = x._test()
    failures, tests = testit
    if tests == 0: continue
    writeln("%s %3d tests, %d failures" % (x.__name__, tests-pt, failures-pf))
    pf, pt = failures, tests

doctest.master.summarize(1)
