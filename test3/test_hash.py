# Hash test for Python 3.2.

import gmpy
import fractions

import sys
try:
    m = sys.hash_info.modulus
except NameError:
    print("new-style hash is not supported")
    sys.exit(0)

for s in [0, m//2, m, m*2, m*m]:
    for i in range(-10,10):
        for k in [-1, 1, 7, 11, -(2**15), 2**16, 2**30, 2**31, 2**32, 2**33, 2**61, -(2**62), 2**63, 2**64]:
            val = k*(s + i)
            assert hash(val) == hash(gmpy.mpz(val))

print("hash tests for integer values passed")

for d in [1, -2, 3, -47, m, m*2]:
    for s in [0, m//2, m, m*2, m*m]:
        for i in range(-10,10):
            for k in [-1, 1, 7, 11, -(2**15), 2**16, 2**30, 2**31, 2**32, 2**33, 2**61, -(2**62), 2**63, 2**64]:
                val = k*(s + i)
                if val:
                    assert hash(fractions.Fraction(d,val)) == hash(gmpy.mpq(d,val)), (d,val,hash(fractions.Fraction(d,val)),hash(gmpy.mpq(d,val)))
                if d:
                    assert hash(fractions.Fraction(val,d)) == hash(gmpy.mpq(val,d)), (val,d,hash(fractions.Fraction(val,d)),hash(gmpy.mpq(val,d)))


print("hash tests for rational values passed")
