import gmpy, time


try: sum
except NameError:
  def sum(x, z):
    for item in x: z += item
    return z

def timedsum(n, zero):
    start=time.clock()
    tot=zero
    for i in range(n):
        tot+=i
    stend=time.clock()
    return type(zero), stend-start, tot

def timedsum1(n, zero):
    start=time.clock()
    tot=sum(range(n), zero)
    stend=time.clock()
    return type(zero), stend-start, tot

def timedmul(n, one):
    start=time.clock()
    tot=one
    for i in range(n):
        tot*=(i+1)
    stend=time.clock()
    return type(one), stend-start, tot

def test(n=100*1000):
    print "Sum of %d items of various types:" % n
    for z in 0L, 0.0, gmpy.mpz(0), gmpy.mpf(0):
        tip, tim, tot = timedsum(n, z)
        print "    %5.3f %.0f %s" % (tim, float(tot), tip)
    print "Sum of %d items of various types w/2.3 sum builtin:" % n
    for z in 0L, 0.0, gmpy.mpz(0), gmpy.mpf(0):
        tip, tim, tot = timedsum1(n, z)
        print "    %5.3f %.0f %s" % (tim, float(tot), tip)
    print "Mul of %d items of various types:" % (n//5)
    for z in 1L, 1.0, gmpy.mpz(1), gmpy.mpf(1):
        tip, tim, tot = timedmul(n//5, z)
        print "    %5.3f %s" % (tim, tip)


if __name__=='__main__':
    test()

