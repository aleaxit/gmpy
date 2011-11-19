# partial unit test for gmpy mpz functionality
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> filter(lambda x: not x.startswith('_'), dir(_g))
['binary', 'bincoef', 'bit_length', 'cdivmod', 'ceil', 'comb', 'denom', 'digits', 'divexact', 'divm', 'f2q', 'fac', 'fbinary', 'fdigits', 'fdivmod', 'fib', 'floor', 'fround', 'fsign', 'fsqrt', 'gcd', 'gcdext', 'get_cache', 'getbit', 'getprec', 'getrprec', 'gmp_limbsize', 'gmp_version', 'hamdist', 'invert', 'is_power', 'is_prime', 'is_square', 'jacobi', 'kronecker', 'lcm', 'legendre', 'license', 'lowbits', 'mpf', 'mpir_version', 'mpq', 'mpz', 'next_prime', 'numdigits', 'numer', 'pi', 'popcount', 'qbinary', 'qdigits', 'qdiv', 'qsign', 'rand', 'reldiff', 'remove', 'root', 'scan0', 'scan1', 'set_cache', 'set_debug', 'set_fcoform', 'set_minprec', 'set_tagoff', 'setbit', 'sign', 'sqrt', 'sqrtrem', 'tdivmod', 'trunc', 'version']
>>> filter(lambda x: not x.startswith('__'), dir(a))
['_copy', 'binary', 'bincoef', 'bit_length', 'comb', 'digits', 'getbit', 'hamdist', 'invert', 'is_power', 'is_prime', 'is_square', 'jacobi', 'kronecker', 'legendre', 'lowbits', 'next_prime', 'numdigits', 'popcount', 'qdiv', 'remove', 'root', 'scan0', 'scan1', 'setbit', 'sign', 'sqrt', 'sqrtrem']
>>>
'''
import warnings
warnings.filterwarnings('ignore', 'setprec')

import gmpy as _g, doctest, sys, operator, gc
__test__={}
a=_g.mpz(123)
b=_g.mpz(456)

# Disable tests since they are not reliable with Python 3.1 but left behind
# in case it is needed in the future.

if sys.platform in ('__DISABLE__linux2', '__DISABLE__darwin'):
  def _memsize():
    """ this function tries to return a measurement of how much memory
        this process is consuming (if it doesn't manage to, it returns 0).
    """
    import os
    try: x = int(os.popen('ps -p %d -o vsz|tail -1' % os.getpid()).read())
    except: x = 0
    return x
else:
  def _memsize():
    return 0

def factorize(x):
    r'''
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

if sys.version[:3] >= '2.5':
  __test__['index']=\
r'''
>>> range(333)[a]
123
>>> range(333)[b]
Traceback (innermost last):
  ...
IndexError: list index out of range
>>> ix = operator.index
>>> ix(_g.mpz(sys.maxint)) == sys.maxint
True
>>> type(ix(_g.mpz(sys.maxint))) is int
True
>>> ix(_g.mpz(sys.maxint+1)) == sys.maxint+1
True
>>> type(ix(_g.mpz(sys.maxint+1))) is long
True
'''

__test__['elemop']=\
r'''
>>> a+b
mpz(579)
>>> a-b
mpz(-333)
>>> a*b
mpz(56088)
>>> a/b
mpz(0)
>>> a%b
mpz(123)
>>> 0%b
mpz(0)
>>> b+a
mpz(579)
>>> b-a
mpz(333)
>>> b*a
mpz(56088)
>>> b%a
mpz(87)
>>> divmod(a,b)
(mpz(0), mpz(123))
>>> divmod(b,a)
(mpz(3), mpz(87))
>>> divmod(0,b)
(mpz(0), mpz(0))
>>> -a
mpz(-123)
>>> a+1
mpz(124)
>>> a+(-1)
mpz(122)
>>> (-1)+a
mpz(122)
>>> 1+a
mpz(124)
>>> a-1
mpz(122)
>>> a-(-1)
mpz(124)
>>> 1-a
mpz(-122)
>>> (-1)-a
mpz(-124)
>>> a+True
mpz(124)
>>> a+False
mpz(123)
>>> a*False
mpz(0)
>>> a//True
mpz(123)
>>> abs(-a)==a
1
>>> pow(a,10)
mpz(792594609605189126649L)
>>> pow(a,7,b)
mpz(99)
>>> _g.sign(b-a)
1
>>> _g.sign(b-b)
0
>>> _g.sign(a-b)
-1
>>> a.sign()
1
>>> (-a).sign()
-1
>>> z=b-b; z.sign()
0
>>> _g.mpz(4)**(0.5)
mpf('2.e0')
>>> import pickle
>>> pickle.loads(pickle.dumps(_g.mpz(12346789)))
mpz(12346789)
>>> s='12345678901234567890123456789'
>>> int(s) == _g.mpz(s)
True
>>> _g.mpz(s) == int(s)
True
>>> del s
'''

from gmpy_truediv import truediv
__test__['newdiv']=\
r'''
>>> a/b
mpz(0)
>>> a//b
mpz(0)
>>> truediv(a,b)
mpf('2.69736842105263157895e-1')
>>> b/a
mpz(3)
>>> b//a
mpz(3)
>>> truediv(b,a)
mpf('3.70731707317073170732e0')
>>>
'''

__test__['divexact']=\
r'''
>>> a=_g.mpz('1234567912345678912345679')
>>> b=_g.mpz('789789789789789789789789')
>>> c=a*b
>>> _g.divexact(c,a)
mpz(789789789789789789789789L)
>>> _g.divexact(10,0)
Traceback (most recent call last):
  ...
ZeroDivisionError: divexact() division by 0
>>>
'''

__test__['divmod']=\
r'''
>>> _g.cdivmod(17,5)
(mpz(4), mpz(-3))
>>> _g.cdivmod(-17,5)
(mpz(-3), mpz(-2))
>>> _g.cdivmod(17,-5)
(mpz(-3), mpz(2))
>>> _g.cdivmod(-17,-5)
(mpz(4), mpz(3))
>>> _g.fdivmod(17,5)
(mpz(3), mpz(2))
>>> _g.fdivmod(-17,5)
(mpz(-4), mpz(3))
>>> _g.fdivmod(17,-5)
(mpz(-4), mpz(-3))
>>> _g.fdivmod(-17,-5)
(mpz(3), mpz(-2))
>>> _g.tdivmod(17,5)
(mpz(3), mpz(2))
>>> _g.tdivmod(-17,5)
(mpz(-3), mpz(-2))
>>> _g.tdivmod(17,-5)
(mpz(-3), mpz(2))
>>> _g.tdivmod(-17,-5)
(mpz(3), mpz(-2))
>>> _g.cdivmod(10,0)
Traceback (most recent call last):
  ...
ZeroDivisionError: cdivmod() division by 0
>>> _g.fdivmod(10,0)
Traceback (most recent call last):
  ...
ZeroDivisionError: fdivmod() division by 0
>>> _g.tdivmod(10,0)
Traceback (most recent call last):
  ...
ZeroDivisionError: tdivmod() division by 0
'''

__test__['cmpr']=\
r'''
>>> c=_g.mpz(a)
>>> c is a
1
>>> c==a
1
>>> c>a
0
>>> c<a
0
>>> d=a._copy()
>>> a is d
0
>>> a == d
1
>>> cmp(a,c)
0
>>> cmp(a,b)
-1
>>> a>b
0
>>> a<b
1
>>> not _g.mpz(0)
1
>>> not a
0
>>> _g.mpz(1) == None
False
>>> _g.mpz(1) == '1'
False
>>> _g.mpz(1) == 'abc'
False
>>> [_g.mpz(23), None].count(None)
1
>>> _g.mpz(3.14)
mpz(3)
>>> _g.mpz(_g.mpq(17,3))
mpz(5)
>>> _g.mpz(23L)
mpz(23)
>>> _g.mpz(-23L)
mpz(-23)
>>> x=1000L*1000*1000*1000*1000*1000*1000
>>> _g.mpz(x)
mpz(1000000000000000000000L)
>>> try: print cmp(_g.mpz(1), _g.mpz(-1))
... except: print 'ouch!'
1
'''

__test__['special'] = \
r'''
>>> a == float('Inf')
False
>>> a != float('Inf')
True
>>> a > float('Inf')
False
>>> a >= float('Inf')
False
>>> a < float('Inf')
True
>>> a <= float('Inf')
True
>>> a == float('-Inf')
False
>>> a != float('-Inf')
True
>>> a > float('-Inf')
True
>>> a >= float('-Inf')
True
>>> a < float('-Inf')
False
>>> a <= float('-Inf')
False
>>> a == float('nan')
False
>>> a != float('nan')
True
>>> a > float('nan')
False
>>> a >= float('nan')
False
>>> a < float('nan')
False
>>> a <= float('nan')
False
>>> float('Inf') == a
False
>>> float('Inf') != a
True
>>> float('Inf') > a
True
>>> float('Inf') >= a
True
>>> float('Inf') < a
False
>>> float('Inf') <= a
False
>>> float('-Inf') == a
False
>>> float('-Inf') != a
True
>>> float('-Inf') > a
False
>>> float('-Inf') >= a
False
>>> float('-Inf') < a
True
>>> float('-Inf') <= a
True
>>> float('nan') == a
False
>>> float('nan') != a
True
>>> float('nan') > a
False
>>> float('nan') >= a
False
>>> float('nan') < a
False
>>> float('nan') <= a
False
>>> a + float('Inf')
inf
>>> float('Inf') + a
inf
>>> a + float('-Inf')
-inf
>>> float('-Inf') + a
-inf
>>> a + float('nan')
nan
>>> float('nan') + a
nan
>>> a - float('Inf')
-inf
>>> float('Inf') - a
inf
>>> a - float('-Inf')
inf
>>> float('-Inf') - a
-inf
>>> a - float('nan')
nan
>>> float('nan') - a
nan
>>> a * float('Inf')
inf
>>> float('Inf') * a
inf
>>> a * float('-Inf')
-inf
>>> float('-Inf') * a
-inf
>>> -a * float('Inf')
-inf
>>> float('Inf') * -a
-inf
>>> -a * float('-Inf')
inf
>>> float('-Inf') * -a
inf
>>> a * float('nan')
nan
>>> float('nan') * a
nan
>>> _g.mpz(0) * float('Inf')
nan
>>> _g.mpz(0) * float('-Inf')
nan
>>> float('Inf') * _g.mpz(0)
nan
>>> float('-Inf') * _g.mpz(0)
nan
>>> a / float('Inf')
mpf('0.e0')
>>> -a / float('Inf')
mpf('0.e0')
>>> float('Inf') / a
inf
>>> float('Inf') / -a
-inf
>>> a / float('-Inf')
mpf('0.e0')
>>> -a / float('-Inf')
mpf('0.e0')
>>> float('-Inf') / a
-inf
>>> float('-Inf') / -a
inf
>>> a / float('nan')
nan
>>> float('nan') / a
nan
>>> float('nan') / _g.mpz(0)
Traceback (most recent call last):
  ...
ZeroDivisionError: mpz division by zero
>>> divmod(a, float('Inf'))
(mpf('0.e0'), mpf('1.23e2'))
>>> divmod(a, float('-Inf'))
(mpf('-1.e0'), -inf)
>>> divmod(-a, float('Inf'))
(mpf('-1.e0'), inf)
>>> divmod(-a, float('-Inf'))
(mpf('0.e0'), mpf('-1.23e2'))
>>> divmod(a, float('nan'))
(nan, nan)
>>> divmod(-a, float('nan'))
(nan, nan)
>>> divmod(_g.mpz(0), float('Inf'))
(mpf('0.e0'), mpf('0.e0'))
>>> divmod(_g.mpz(0), float('-Inf'))
(mpf('0.e0'), mpf('0.e0'))
>>> divmod(_g.mpz(0), float('nan'))
(nan, nan)
>>> divmod(float('Inf'), a)
(nan, nan)
>>> divmod(float('-Inf'), a)
(nan, nan)
>>> divmod(float('Inf'), -a)
(nan, nan)
>>> divmod(float('-Inf'), -a)
(nan, nan)
>>> divmod(float('nan'), a)
(nan, nan)
>>> divmod(float('nan'), -a)
(nan, nan)
>>> divmod(float('Inf'), _g.mpz(0))
Traceback (most recent call last):
  ...
ZeroDivisionError: mpz modulo by zero
>>> divmod(float('-Inf'), _g.mpz(0))
Traceback (most recent call last):
  ...
ZeroDivisionError: mpz modulo by zero
>>> divmod(float('nan'), _g.mpz(0))
Traceback (most recent call last):
  ...
ZeroDivisionError: mpz modulo by zero
'''

__test__['bitops']=\
r'''
>>> ~a
mpz(-124)
>>> a&b
mpz(72)
>>> a|b
mpz(507)
>>> a^b
mpz(435)
>>> a<<1
mpz(246)
>>> a>>1
mpz(61)
>>> a<<-1
Traceback (innermost last):
  ...
ValueError: negative shift count
>>> a>>-2
Traceback (innermost last):
  ...
ValueError: negative shift count
>>> a<<0
mpz(123)
>>> a>>0
mpz(123)
>>> a.popcount()
6
>>> _g.popcount(b)
4
>>> _g.popcount(-7)
-1
>>> _g.popcount(0)
0
>>> a.hamdist(b)
6
>>> _g.hamdist(3)
Traceback (innermost last):
  ...
TypeError: hamdist() expects 'mpz','mpz' arguments
>>> a.hamdist()
Traceback (innermost last):
  ...
TypeError: hamdist() expects 'mpz','mpz' arguments
>>> a.hamdist(3, 4)
Traceback (innermost last):
  ...
TypeError: hamdist() expects 'mpz','mpz' arguments
>>> a.lowbits(5)
mpz(27)
>>> b.lowbits(5)
mpz(8)
>>> b.lowbits(5)==(b%32)
1
>>> a.lowbits(5)==(a%32)
1
>>> a.setbit(20)
mpz(1048699)
>>> a.setbit(0,0)
mpz(122)
>>> for i in range(8):
...     print a.getbit(i),
...     if i==7: print
...
1 1 0 1 1 1 1 0
>>> for i in range(10):
...     print b.getbit(i),
...     if i==9: print
...
0 0 0 1 0 0 1 1 1 0
>>> [a.scan0(j) for j in range(33)]
[2, 2, 2, 7, 7, 7, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
>>> [a.scan1(j) for j in range(10)]
[0, 1, 3, 3, 4, 5, 6, None, None, None]
>>> n=_g.mpz(-(7+6*16+5*256+7*4092))
>>> [n.scan0(j) for j in range(18)]
[1, 1, 3, 3, 6, 6, 6, 8, 8, 10, 10, 12, 12, 13, 14, -1, None, None]
>>> [n.scan1(j) for j in range(33)]
[0, 2, 2, 4, 4, 5, 7, 7, 9, 9, 11, 11, 15, 15, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
>>> _g.mpz(0).bit_length()
0
>>> _g.mpz(12345).bit_length()
14
'''

__test__['format']=\
r'''
>>> str(a)
'123'
>>> repr(a)
'mpz(123)'
>>> hex(a)
'0x7b'
>>> oct(a)
'0173'
>>> _g.mpz('123')
mpz(123)
>>> _g.mpz('1001001011',2)
mpz(587)
>>> _g.mpz('1001001011',2).digits(2)
'1001001011'
>>> for i in range(2,63):
...     print a.digits(i),
...     if i%6==0: print
...
1111011 11120 1323 443 323
234 173 146 123 102 a3
96 8b 83 7b 74 6f
69 63 5i 5d 58 53
4n 4j 4f 4b 47 43
3u 3r 3o 3l 3i 3f
3C 39 36 33 30 2d
2b 2Z 2X 2V 2T 2R
2P 2N 2L 2J 2H 2F
2D 2B 29 27 25 23
21 1z
>>> print a.digits(63)
Traceback (innermost last):
  File "<stdin>", line 1, in ?
    print a.digits(37)
ValueError: base must be either 0 or in the interval 2 ... 62
>>> _g.set_tagoff(0)
1
>>> a
gmpy.mpz(123)
>>> _g.set_tagoff(1)
0
>>> _g.mpz('43')
mpz(43)
>>> _g.mpz('043')
mpz(43)
>>> _g.mpz('43',0)
mpz(43)
>>> _g.mpz('043',0)
mpz(35)
>>> _g.mpz('0x43',0)
mpz(67)
>>> _g.mpz('0x43')
Traceback (innermost last):
  File "<pyshell#181>", line 1, in ?
    _g.mpz('0x43')
ValueError: invalid digits
>>>
'''

__test__['binio']=\
r'''
>>> ba=a.binary()
>>> ba
'{'
>>> _g.mpz(ba,256)
mpz(123)
>>> _g.mpz(ba,256)==a
1
>>> _g.binary(123)
'{'
>>> z=_g.mpz('melancholy',256)
>>> z
mpz(573406620562849222387053L)
>>> long(z)
573406620562849222387053L
>>> divmod(z,a)
(mpz(4661842443600400182008L), mpz(69))
>>> for i in range(2,63):
...    print i,z.numdigits(i),
...    if i%6==0: print
...
2 79 3 50 4 40 5 35 6 31
7 29 8 27 9 25 10 24 11 23 12 23
13 22 14 21 15 21 16 20 17 20 18 19
19 19 20 19 21 18 22 18 23 18 24 18
25 18 26 17 27 17 28 17 29 17 30 17
31 16 32 16 33 16 34 16 35 16 36 16
37 16 38 16 39 15 40 15 41 15 42 15
43 15 44 15 45 15 46 15 47 15 48 15
49 15 50 14 51 14 52 14 53 14 54 14
55 14 56 14 57 14 58 14 59 14 60 14
61 14 62 14
>>> _g.numdigits(23)
2
>>> _g.numdigits(23,2)
5
>>> _g.numdigits(23,99)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be either 0 or in the interval 2 ... 62
>>> hash(a)
123
>>> hash(b)
456
>>> hash(z) == hash(long(z))
True
>>> _g.mpz(_g.binary(-123),256)
mpz(-123)
>>> long(_g.mpz(-3))
-3L
>>>
'''

__test__['number']=\
r'''
>>> print a.sqrt(), b.sqrt()
11 21
>>> print a.sqrtrem(), b.sqrtrem()
(mpz(11), mpz(2)) (mpz(21), mpz(15))
>>> for i in range(5):
...    print a.root(i+1),b.root(i+1)
...
(mpz(123), 1) (mpz(456), 1)
(mpz(11), 0) (mpz(21), 0)
(mpz(4), 0) (mpz(7), 0)
(mpz(3), 0) (mpz(4), 0)
(mpz(2), 0) (mpz(3), 0)
>>> a.is_square()
0
>>> a.is_power()
0
>>> _g.is_square(99*99)
1
>>> _g.is_square(99*99*99)
0
>>> _g.is_square(0)
1
>>> _g.is_square(-1)
0
>>> _g.is_power(99*99*99)
1
>>> _g.gcd(a,b)
mpz(3)
>>> temp=_g.gcdext(a,b)
>>> temp[0]==a*temp[1]+b*temp[2]
True
>>> _g.lcm(a,b)
mpz(18696)
>>> _g.fac(7)
mpz(5040)
>>> _g.fib(17)
mpz(1597)
>>> for i in range(10):
...     print _g.bincoef(10,i)
...
1
10
45
120
210
252
210
120
45
10
>>> _g.divm(b,a,20)
mpz(12)
>>> _g.divm(a,b,100)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.divm(a,b,100)
ZeroDivisionError: not invertible
>>> _g.divm(6,12,14)
mpz(4)
>>> _g.divm(0,1,2)
mpz(0)
>>> # guard against regression of an ancient gmpy bug: divm w/non-coprime parms
>>> _g.divm(4,8,20)
mpz(3)
>>> _g.divm(4,8,20)
mpz(3)
>>> _g.mpz(20)
mpz(20)
>>> _g.mpz(8)
mpz(8)
>>> _g.mpz(4)
mpz(4)
>>> # guard against regression of a memory leak in divm
>>> __ = gc.collect()
>>> _siz = 87654
>>> _siz = _memsize()
>>> for x in xrange(45678):
...     _xx=_g.divm(b,a,20)
>>> del _xx
>>> __ = gc.collect()
>>> (_memsize()-_siz) <= 16
True
>>> a.invert(100)
mpz(87)
>>> b.invert(100)
mpz(0)
>>> _g.invert(3)
Traceback (innermost last):
  ...
TypeError: invert() expects 'mpz','mpz' arguments
>>> a.invert()
Traceback (innermost last):
  ...
TypeError: invert() expects 'mpz','mpz' arguments
>>> a.invert(3, 4)
Traceback (innermost last):
  ...
TypeError: invert() expects 'mpz','mpz' arguments
>>> _g.comb(3,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: binomial coefficient with negative k
>>> _g.sqrt(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: sqrt of negative number
>>> _g.sqrtrem(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: sqrt of negative number
>>> _g.remove(3,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: factor must be > 1
>>> _g.remove(3)
Traceback (innermost last):
  ...
TypeError: remove() expects 'mpz','mpz' arguments
>>> a.remove()
Traceback (innermost last):
  ...
TypeError: remove() expects 'mpz','mpz' arguments
>>> a.remove(3, 4)
Traceback (innermost last):
  ...
TypeError: remove() expects 'mpz','mpz' arguments
>>> _g.is_prime(3,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: repetition count for is_prime must be positive
>>> _g.jacobi(10,3)
1
>>> _g.jacobi(10,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: jacobi's y must be odd prime > 0
>>> _g.jacobi(3)
Traceback (innermost last):
  ...
TypeError: jacobi() expects 'mpz','mpz' arguments
>>> a.jacobi()
Traceback (innermost last):
  ...
TypeError: jacobi() expects 'mpz','mpz' arguments
>>> a.jacobi(3, 4)
Traceback (innermost last):
  ...
TypeError: jacobi() expects 'mpz','mpz' arguments
>>> _g.legendre(10,3)
1
>>> _g.legendre(10,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: legendre's y must be odd and > 0
>>> _g.legendre(3)
Traceback (innermost last):
  ...
TypeError: legendre() expects 'mpz','mpz' arguments
>>> a.legendre()
Traceback (innermost last):
  ...
TypeError: legendre() expects 'mpz','mpz' arguments
>>> a.legendre(3, 4)
Traceback (innermost last):
  ...
TypeError: legendre() expects 'mpz','mpz' arguments
>>> _g.kronecker(10,3)
1
>>> _g.kronecker(10,-3)
1
>>> _g.kronecker(3)
Traceback (innermost last):
  ...
TypeError: kronecker() expects 'mpz','mpz' arguments
>>> a.kronecker()
Traceback (innermost last):
  ...
TypeError: kronecker() expects 'mpz','mpz' arguments
>>> a.kronecker(3, 4)
Traceback (innermost last):
  ...
TypeError: kronecker() expects 'mpz','mpz' arguments
>>> a=10L**20
>>> b=a+39
>>> _g.jacobi(a,b)
1
>>> _g.legendre(a,b)
1
>>> _g.kronecker(a,b)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: Either arg in Kronecker must fit in an int
>>> f=_g.mpf(3.3)
>>> f.setprec(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: n must be >=0
>>> _g.rand('init',-1)
>>> _g.rand('init',-7)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: size must be in 1..128
>>> _g.rand('init',200)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: size must be in 1..128
>>> _g.rand('qual')
32
>>> _g.rand('floa',-7)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: 'floa' needs arg>=0
'''

def _test(chat=None):
    try:
        float('nan')
    except ValueError:
        if "special" in __test__:
            del(__test__["special"])

    if chat:
        print "Unit tests for gmpy 1.15 (mpz functionality)"
        print "    running on Python %s" % sys.version
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
    _g.set_cache(0,128)

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
        print "Overall results for mpz:"
    return doctest.master.summarize(chat)


if __name__=='__main__':
    _test(1)

