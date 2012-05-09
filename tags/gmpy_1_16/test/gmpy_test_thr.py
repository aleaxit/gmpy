# partial unit test for gmpy threaded mpz functionality
# relies on Tim Peters' "doctest.py" test-driver


import gmpy as _g, doctest, sys, operator, gc, Queue, threading

__test__={}
def _tf(N=2, _K=1234**5678):
    """Takes about 100ms on a first-generation Macbook Pro"""
    for i in range(N): assert (_g.mpz(1234)**5678)==_K
a=_g.mpz(123)
b=_g.mpz(456)
c=_g.mpz(123456789123456789)

def factorize(x=c):
    r'''
    (Takes about 25ms, on c, on a first-generation Macbook Pro)
    >>> factorize(a)
    [3, 41]
    >>> factorize(b)
    [2, 2, 2, 3, 19]
    >>>
    '''
    import gmpy as _g
    savex=x
    prime=2
    x=_g.mpz(x)
    factors=[]
    while x>=prime:
        newx,mult=x.remove(prime)
        if mult:
            factors.extend([int(prime)]*mult)
            x=newx
        prime=_g.next_prime(prime)
    for factor in factors: assert _g.is_prime(factor)
    from operator import mul
    assert reduce(mul, factors)==savex
    return factors

def elemop(N=1000):
    r'''
    (Takes about 40ms on a first-generation Macbook Pro)
    '''
    for i in range(N):
        assert a+b == 579
        assert a-b == -333
        assert b*a == a*b == 56088
        assert b%a == 87
        assert divmod(a, b) == (0, 123)
        assert divmod(b, a) == (3, 87)
        assert -a == -123
        assert pow(a, 10) == 792594609605189126649
        assert pow(a, 7, b) == 99
        assert cmp(a, b) == -1
        assert '7' in str(c)
        assert '0' not in str(c)
        assert a.sqrt() == 11
        assert _g.lcm(a, b) == 18696
        assert _g.fac(7) == 5040
        assert _g.fib(17) == 1597
        assert _g.divm(b, a, 20) == 12
        assert _g.divm(4, 8, 20) == 3
        assert _g.divm(4, 8, 20) == 3
        assert _g.mpz(20) == 20
        assert _g.mpz(8) == 8
        assert _g.mpz(4) == 4
        assert a.invert(100) == 87

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy 1.15 (threading)"
        print "    running on Python", sys.version
        print
        if _g.gmp_version():
            print "Testing gmpy %s (GMP %s) with default caching (%s, %s)" % (
                (_g.version(), _g.gmp_version(), _g.get_cache()[0],
                _g.get_cache()[1]))
        else:
            print "Testing gmpy %s (MPIR %s) with default caching (%s, %s)" % (
                (_g.version(), _g.mpir_version(), _g.get_cache()[0],
                _g.get_cache()[1]))

    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat: print "Repeating tests, with caching disabled"
    _g.set_zcache(0)

    sav = sys.stdout
    class _Dummy:
        def write(self,*whatever):
            pass
    try:
        sys.stdout = _Dummy()
        doctest.testmod(thismod, report=0)
    finally:
        sys.stdout = sav

    if chat:
        print
        print "Overall results for thr:"
    return doctest.master.summarize(chat)

class DoOne(threading.Thread):
    def __init__(self, q):
        threading.Thread.__init__(self)
        self.q = q
    def run(self):
        while True:
            task = self.q.get()
            if task is None: break
            task()

def _test_thr(Ntasks=5, Nthreads=1):
    q = Queue.Queue()
    funcs = (_tf, 1), (factorize, 4), (elemop, 2)
    for i in range(Ntasks):
        for f, n in funcs:
            for x in range(n):
                q.put(f)
    for i in range(Nthreads):
        q.put(None)
    thrs = [DoOne(q) for i in range(Nthreads)]
    for t in thrs: t.start()
    for t in thrs: t.join()

if __name__=='__main__':
    _test(1)

