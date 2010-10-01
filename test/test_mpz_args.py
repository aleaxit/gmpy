# Test a wide variety of input values to the commonly used mpz operations.
# This test should be run whenever optimizations are made to the handling of
# arguments.

import sys
import gmpy

if sys.version.startswith('3'):
    intTypes = (int,)
else:
    intTypes = (int, long)

def writeln(s):
    sys.stdout.write(s+'\n')

valueList = [0, 1, 2, 3, 4, 5]

for power in (15, 16, 30, 32, 45, 48, 60, 64, 75, 90, 96, 105, 120, 128):
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
    mpzValues.append(gmpy.mpz(i))
    mpzValues.append(-gmpy.mpz(i))

testValues.extend(mpzValues)

for i in testValues:
    for z in mpzValues:
        # Test all permutations of addition
        assert int(i)+int(z) == i+z, (repr(i),repr(z))
        assert int(z)+int(i) == z+i, (repr(i),repr(z))

        # Test all permutations of subtraction
        assert int(i)-int(z) == i-z, (repr(i),repr(z))
        assert int(z)-int(i) == z-i, (repr(i),repr(z))

        # Test all permutations of multiplication
        assert int(i)*int(z) == i*z, (repr(i),repr(z))
        assert int(z)*int(i) == z*i, (repr(i),repr(z))

        # Test all permutations of division
        if z!=0:
            temp = int(i)//int(z)
            assert int(i)//int(z) == i//z, (repr(i),repr(z))
            assert int(i)%int(z) == i%z, (repr(i),repr(z))
            assert divmod(int(i),int(z)) == divmod(i,z), (repr(i),repr(z))

        if i!=0:
            temp = int(z)//int(i)
            assert int(z)//int(i) == z//i, (repr(i),repr(z))
            assert int(z)%int(i) == z%i, (repr(i),repr(z))
            assert divmod(int(z),int(i)) == divmod(z,i), (repr(i),repr(z))
