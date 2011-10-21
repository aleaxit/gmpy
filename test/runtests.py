from __future__ import print_function

import sys
import gmpy2

print()
print("Unit tests for gmpy2 {0} with Python {1}".format(gmpy2.version(), sys.version.split()[0]))
print("  Mutliple-precision library:   {0}".format(gmpy2.mp_version()))
print("  Floating-point library:       {0}".format(gmpy2.mpfr_version()))
print("  Complex library:              {0}".format(gmpy2.mpc_version()))
print("  Caching Values: (Number)      {0}".format(gmpy2.get_cache()[0]))
print("  Caching Values: (Size, limbs) {0}".format(gmpy2.get_cache()[1]))
print()
print("No tests found!")
