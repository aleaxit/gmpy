from __future__ import print_function
import gmpy2

def test_pack_unpack(repeat):
    """Test gmpy2.pack and gmpy2.unpack."""
    r = gmpy2.random_state(42)
    for counter in range(repeat):
        for t in (10, 1000, 2000, 10000, 100000):
            v = gmpy2.mpz_rrandomb(r, t)
            for b in range(1, max(1001,t)):
                temp = gmpy2.unpack(v, b)
                u = gmpy2.pack(temp, b)
                assert u == v


if __name__ == "__main__":
    import sys
    try:
        repeat = abs(int(sys.argv[1]))
    except:
        repeat = 1
    print("Testing pack/unpack for a large number of values.")
    print("This test may take a few minutes.")
    test_pack_unpack(repeat)
    print("Test successful.")
