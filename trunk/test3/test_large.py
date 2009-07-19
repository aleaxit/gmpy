from gmpy import *
from math import log
from time import time

# This test is designed to detect issues when allocating memory for large
# numbers. If it crashes and you need to work with very large numbers,
# you will need to compile GMP from scratch and try a different memory
# allocation option.



def pi(N):
    print("Computing pi to %s decimal places." % N)
    start = time()
    N = int(round(log(10,2)*N))
    sq2 = fsqrt(mpf(2, N))
    a = mpz(6) - 4*sq2
    y = sq2-1
    for k in range(0, 10000):
        xx = fsqrt(fsqrt(1-y**4))
        y = (1-xx)/(1+xx)
        anew = a*(1+y)**4 - 2**(2*k+3)*y*(1+y+y**2)
        if anew == a:
            break
        a = anew
    print("Computation took %5.2f seconds." % (time() - start))
    return 1/a

if __name__ == '__main__':
    print("Testing operations with large numbers.")
    pi(1000000)
