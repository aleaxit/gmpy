from __future__ import print_function
import gmpy2

def test_pack_unpack(bits = 200, chunk = 500, terms = 50):
    """Test gmpy2.pack and gmpy2.unpack."""
    for t in range(2, terms):
        print(t)
        for b in range(1, bits):
            # Test with all bits set to 1.
            v = [ 2**b - 1 ] * t
            for c in range(b, chunk):
                temp = gmpy2.pack(v, c)
                u = gmpy2.unpack(temp, c)
                assert u == v, (v, temp, u, (t, b, c))

def main():
    test_pack_unpack()

if __name__ == "__main__":
    main()
