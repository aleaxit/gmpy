import gmpy, time

def timedfib(n, zero):
    start=time.clock()
    a=zero; b=a+1
    for i in range(n):
        a,b=b,a+b
    stend=time.clock()
    return type(zero), stend-start, a

def timedfibsp(n, one):
    start=time.clock()
    result=gmpy.fib(n)
    stend=time.clock()
    return type(one), stend-start, result

def test(n=100*1000):
    print "%dth Fibonacci number of various types:" % n
    for z in 0L, gmpy.mpz(0), gmpy.mpf(0):
        tip, tim, tot = timedfib(n, z)
        print "    %5.3f %s %s" % (tim, gmpy.fdigits(tot,10,6), tip)
    tip, tim, tot = timedfibsp(n, 1)
    print "    %5.3f %s %s" % (tim, gmpy.fdigits(tot,10,6), "gmpy.fib")


if __name__=='__main__':
    test()

