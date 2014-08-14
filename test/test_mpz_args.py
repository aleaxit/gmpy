# Test a wide variety of input values to the commonly used mpz operations.
# This test should be run whenever optimizations are made to the handling of
# arguments.

import sys
import gmpy2

if sys.version.startswith('3'):
    intTypes = (int,)
else:
    intTypes = (int, long)


valueList = [0, 1, 2, 3, 4, 5]

for power in (14, 15, 16, 29, 30, 31, 32, 45, 48, 60, 64, 75, 90, 96, 105, 120, 128):
    for i in (-2, -1, 0, 1, 2):
        valueList.append(2**power + i)

valueList.append('123456789012345678901234567890')
valueList.append('10000000000000000000000000000000000000000000000000000000000000000')

testValues = []
mpzValues = []
for i in valueList:
    for t in intTypes:
        testValues.append(t(i))
        testValues.append(-t(i))
    mpzValues.append(gmpy2.mpz(i))
    mpzValues.append(-gmpy2.mpz(i))

testValues.extend(mpzValues)

def test():
    for i in testValues:
        for z in mpzValues:
            # Test all permutations of addition
            assert int(i) + int(z) == i + z, (repr(i), repr(z))
            assert int(z) + int(i) == z + i, (repr(i), repr(z))

            # Test all permutations of subtraction
            assert int(i) - int(z) == i - z, (repr(i), repr(z))
            assert int(z) - int(i) == z - i, (repr(i), repr(z))

            # Test all permutations of multiplication
            assert int(i) * int(z) == i * z, (repr(i), repr(z))
            assert int(z) * int(i) == z * i, (repr(i), repr(z))

            # Test all permutations of floor division
            if z != 0:
                assert int(i) // int(z) == i // z, (repr(i), repr(z))
                assert int(i) % int(z) == i % z, (repr(i), repr(z))
                assert divmod(int(i), int(z)) == divmod(i, z), (repr(i), repr(z))

            if i!=0:
                assert int(z) // int(i) == z // i, (repr(i), repr(z))
                assert int(z) % int(i) == z % i, (repr(i), repr(z))
                assert divmod(int(z), int(i)) == divmod(z,i), (repr(i), repr(z))
    return True
            
if __name__ == "__main__":
    print("Testing combinations of mpz and integer operations.")
    print("This test may take a few minutes.")
    test()
    print("Test successful.")
