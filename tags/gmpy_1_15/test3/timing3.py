import gmpy as _g
import time

print("Typical expected results would be:","""
D:\PySym>python timing.py
Factorial of 10000 took 0.0619989238859 (35660 digits)
Fibonacci of 10000 took 0.000744228458022 (2090 digits)
Factorial of 100000 took 4.44311764676 (456574 digits)
Fibonacci of 100000 took 0.022344453738 (20899 digits)
Factorial of 1000000 took 152.151135367 (5565709 digits)
Fibonacci of 1000000 took 0.670207059778 (208988 digits)
""")

print("Actual timings and results...:")
for i in (10000,100000,1000000):
    start=time.time()
    x=_g.fac(i)
    stend=time.time()
    print("Factorial of %d took %s (%d digits)" % (
        i, stend-start, x.numdigits()))

    start=time.time()
    x=_g.fib(i)
    stend=time.time()
    print("Fibonacci of %d took %s (%d digits)" % (
        i, stend-start, x.numdigits()))

