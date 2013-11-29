from __future__ import print_function
import gmpy2
import fractions
import decimal

MPZ = gmpy2.mpz
MPQ = gmpy2.mpq
MPF = gmpy2.mpf
FR = fractions.Fraction
DC = decimal.Decimal

int_vals = []
int_vals.extend([1, 2, 3, -10, 256, 23**65])
int_vals.extend([MPZ(1), MPZ(2), MPZ(3), MPZ(-10), MPZ(23**65)])

float_vals = []
float_vals.extend([1.23, -3.124159, 0.0234, float("nan")])
float_vals.extend([MPF(1.23), MPF("-3.124159"), MPF("0.0234")])

frac_vals = []
frac_vals.extend([FR(1,2), FR(-7,2345)])
frac_vals.extend([MPQ(1,2), MPQ(-7,2345)])

all_vals = int_vals + float_vals + frac_vals

def test_leaks1(bits = 80, chunk = 150, terms = 20):
    """Test gmpy2.pack and gmpy2.unpack."""

    for t in range(2, terms):
        for b in range(1, bits):
            # Test with all bits set to 1.
            v = [ 2**b - 1 ] * t
            for c in range(b, chunk):
                temp = gmpy2.pack(v, c)
                u = gmpy2.unpack(temp, c)
                assert u == v, (v, temp, u, (t, b, c))

def test_leaks2():
    """Test binary operations."""

    def test_binary(a, b):
        t = a + b
        t = a - b
        t = a * b
        try:
            t = a / b
        except:
            pass
        try:
            t = a // b
        except:
            pass
        try:
            t = divmod(a, b)
        except:
            pass
        t = gmpy2.add(a, b)
        t = gmpy2.sub(a, b)
        t = gmpy2.mul(a, b)
        t = gmpy2.div(a, b)
        t = gmpy2.agm(a, b)
        t = gmpy2.atan2(a, b)

    for x in all_vals:
        for y in all_vals:
            test_binary(x, y)
            test_binary(y, x)


def main(count = 100):
    print("Testing for memory leaks by repeating a set of calculations.")
    print("This test may take a few minutes.")
    for i in range(count):
        print("Pass: {0}".format(i))
        test_leaks2()
    print("Test successful.")

if __name__ == "__main__":
    main()
