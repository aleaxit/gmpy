# partial unit test for gmpy2 xmpz functionality
# relies on Tim Peters' "doctest.py" test-driver

import gmpy2 as _g, doctest, sys, operator, gc
__test__={}
a=_g.xmpz(123)
b=_g.xmpz(456)
aa=_g.mpz(123)
bb=_g.mpz(456)

__test__['index']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> range(333)[a]
123
>>> range(333)[b]
Traceback (innermost last):
  ...
IndexError: range object index out of range
'''

__test__['elemop']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> a+b
mpz(579)
>>> a-b
mpz(-333)
>>> a*b
mpz(56088)
>>> a//b
mpz(0)
>>> a/b
mpf('2.69736842105263157895e-1')
>>> b//a
mpz(3)
>>> b/a
mpf('3.70731707317073170732e0')
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
>>> a
xmpz(-123)
>>> a=_g.xmpz(123)
>>> a+1
mpz(124)
>>> a+=1; a
xmpz(124)
>>> a=_g.xmpz(123)
>>> a+(-1)
mpz(122)
>>> a-=1; a
xmpz(122)
>>> a=_g.xmpz(123)
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
>>> -a
>>> abs(a)
>>> a
xmpz(123)
>>> pow(a,10)
mpz(792594609605189126649)
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
>>> -a
>>> a.sign()
-1
>>> a=_g.xmpz(123)
>>> z=b-b; z.sign()
0
>>> s='12345678901234567890123456789'
>>> int(s) == _g.xmpz(s)
True
>>> _g.xmpz(s) == int(s)
True
>>> del s
>>> _g.is_even(a)
False
>>> _g.is_odd(a)
True
>>> _g.is_even(b)
True
>>> _g.is_odd(b)
False
>>> a.is_odd()
True
>>> a.is_even()
False
>>> _g.is_even(2)
True
'''

__test__['special'] = \
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
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
mpf('inf')
>>> float('Inf') + a
mpf('inf')
>>> a + float('-Inf')
mpf('-inf')
>>> float('-Inf') + a
mpf('-inf')
>>> a + float('nan')
mpf('nan')
>>> float('nan') + a
mpf('nan')
>>> a - float('Inf')
mpf('-inf')
>>> float('Inf') - a
mpf('inf')
>>> a - float('-Inf')
mpf('inf')
>>> float('-Inf') - a
mpf('-inf')
>>> a - float('nan')
mpf('nan')
>>> float('nan') - a
mpf('nan')
>>> a * float('Inf')
mpf('inf')
>>> float('Inf') * a
mpf('inf')
>>> a * float('-Inf')
mpf('-inf')
>>> float('-Inf') * a
mpf('-inf')
>>> -a
>>> a * float('Inf')
mpf('-inf')
>>> float('Inf') * a
mpf('-inf')
>>> a * float('-Inf')
mpf('inf')
>>> float('-Inf') * a
mpf('inf')
>>> -a
>>> a * float('nan')
mpf('nan')
>>> float('nan') * a
mpf('nan')
>>> _g.xmpz(0) * float('Inf')
mpf('nan')
>>> _g.xmpz(0) * float('-Inf')
mpf('nan')
>>> float('Inf') * _g.xmpz(0)
mpf('nan')
>>> float('-Inf') * _g.xmpz(0)
mpf('nan')
>>> a
2
>>> a / float('Inf')
mpf('0.0e0')
>>> float('Inf') / a
mpf('inf')
>>> a / float('-Inf')
mpf('-0.0e0')
>>> float('-Inf') / a
mpf('-inf')
>>> a / float('nan')
mpf('nan')
>>> float('nan') / a
mpf('nan')
>>> -a
>>> a / float('Inf')
mpf('0.0e0')
>>> float('Inf') / a
mpf('-inf')
>>> a / float('-Inf')
mpf('0.0e0')
>>> float('-Inf') / a
mpf('inf')
>>> -a
'''

__test__['divexact']=\
r'''
>>> a=_g.xmpz('1234567912345678912345679')
>>> b=_g.xmpz('789789789789789789789789')
>>> c=a*b
>>> _g.divexact(c,a)
mpz(789789789789789789789789)
>>>
'''

__test__['divmod']=\
r'''
>>> _g.set_prefer_mutable(1)
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
>>> _g.set_prefer_mutable(0)
'''

__test__['cmpr']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> c=_g.xmpz(a)
>>> c is a
0
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
>>> a>b
0
>>> a<b
1
>>> not _g.xmpz(0)
1
>>> not a
0
>>> _g.xmpz(1) == None
False
>>> _g.xmpz(1) == '1'
False
>>> _g.xmpz(1) == 'abc'
False
>>> [_g.xmpz(23), None].count(None)
1
>>> _g.xmpz(3.14)
xmpz(3)
>>> _g.xmpz(_g.mpq(17,3))
xmpz(5)
>>> _g.xmpz(23)
xmpz(23)
>>> _g.xmpz(-23)
xmpz(-23)
>>> x=1000*1000*1000*1000*1000*1000*1000
>>> _g.xmpz(x)
xmpz(1000000000000000000000)
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
'''

__test__['bitops']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> ~a
>>> a
xmpz(-124)
>>> a=_g.xmpz(123)
>>> a&b
xmpz(72)
>>> a|b
xmpz(507)
>>> a^b
xmpz(435)
>>> a<<1
xmpz(246)
>>> a>>1
xmpz(61)
>>> a<<-1
Traceback (innermost last):
  ...
ValueError: negative shift count
>>> a>>-2
Traceback (innermost last):
  ...
ValueError: negative shift count
>>> a<<0
xmpz(123)
>>> a>>0
xmpz(123)
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
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.hamdist()
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.hamdist(3, 4)
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.bit_set(20)
>>> a
xmpz(1048699)
>>> a=_g.xmpz(123)
>>> a.bit_clear(0)
>>> a
xmpz(122)
>>> a=_g.xmpz(123)
>>> for i in range(8):
...     print(a.bit_test(i))
...
True
True
False
True
True
True
True
False
>>> for i in range(10):
...     print(b.bit_test(i))
...
False
False
False
True
False
False
True
True
True
False
>>> [a.bit_scan0(j) for j in range(33)]
[2, 2, 2, 7, 7, 7, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
>>> [a.bit_scan1(j) for j in range(10)]
[0, 1, 3, 3, 4, 5, 6, None, None, None]
>>> n=_g.xmpz(-(7+6*16+5*256+7*4092))
>>> [n.bit_scan0(j) for j in range(18)]
[1, 1, 3, 3, 6, 6, 6, 8, 8, 10, 10, 12, 12, 13, 14, -1, None, None]
>>> [n.bit_scan1(j) for j in range(33)]
[0, 2, 2, 4, 4, 5, 7, 7, 9, 9, 11, 11, 15, 15, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
>>> _g.xmpz(0).bit_length()
0
>>> _g.xmpz(12345).bit_length()
14
>>> _old=_g.get_prefer_mutable()
>>> _g.set_prefer_mutable(1)
>>> _g.bit_mask(9)
xmpz(511)
>>> _g.set_prefer_mutable(_old)
>>> del(_old)
'''

__test__['format']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> str(a)
'123'
>>> repr(a)
'xmpz(123)'
>>> hex(a)
'0x7b'
>>> oct(a)
'0o173'
>>> _g.xmpz('123')
xmpz(123)
>>> _g.xmpz('1001001011',2)
xmpz(587)
>>> _g.xmpz('0b1001001011')
xmpz(587)
>>> _g.xmpz('1001001011',2).digits(2)
'0b1001001011'
>>> for i in range(2,37):
...     print(a.digits(i))
...
0b1111011
11120
1323
443
323
234
0o173
146
123
102
a3
96
8b
83
0x7b
74
6f
69
63
5i
5d
58
53
4n
4j
4f
4b
47
43
3u
3r
3o
3l
3i
3f
>>> print(a.digits(37))
Traceback (innermost last):
  File "<stdin>", line 1, in ?
    print(a.digits(37))
ValueError: base must be either 0 or in the interval 2 ... 36
>>> _g.xmpz('43')
xmpz(43)
>>> _g.xmpz('043')
xmpz(43)
>>> _g.xmpz('43',10)
xmpz(43)
>>> _g.xmpz('0o43',0)
xmpz(35)
>>> _g.xmpz('0x43')
xmpz(67)
>>> _g.xmpz('0x43',10)
Traceback (innermost last):
  File "<pyshell#181>", line 1, in ?
    _g.xmpz('0x43')
ValueError: invalid digits
>>>
'''

__test__['binio']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> ba=a.binary()
>>> ba
b'{'
>>> _g.xmpz(ba,256)
xmpz(123)
>>> _g.xmpz(ba,256)==a
1
>>> _g.binary(_g.xmpz(123))
b'{'
>>> z=_g.xmpz('melancholy',256)
>>> z
xmpz(573406620562849222387053)
>>> int(z)
573406620562849222387053
>>> divmod(z,a)
(xmpz(4661842443600400182008), xmpz(69))
>>> for i in range(2,37):
...    print(i,z.numdigits(i))
...
2 79
3 50
4 40
5 35
6 31
7 29
8 27
9 25
10 24
11 23
12 23
13 22
14 21
15 21
16 20
17 20
18 19
19 19
20 19
21 18
22 18
23 18
24 18
25 18
26 17
27 17
28 17
29 17
30 17
31 16
32 16
33 16
34 16
35 16
36 16
>>> _g.numdigits(23)
2
>>> _g.numdigits(23,2)
5
>>> _g.numdigits(23,99)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be either 0 or in the interval 2 ... 36
>>> _g.xmpz(_g.binary(_g.xmpz(-123)),256)
xmpz(-123)
>>> int(_g.xmpz(-3))
-3
'''

__test__['number']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> a.rootrem(2)
(xmpz(11), xmpz(2))
>>> a.rootrem(3)
(xmpz(4), xmpz(59))
>>> _g.rootrem(a*a)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: rootrem() requires 'mpz','int' arguments
>>> _g.rootrem(a*a,2)
(xmpz(123), xmpz(0))
>>> a.sqrt()
>>> print(a)
11
>>> b.sqrt()
>>> print(b)
21
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
>>> print(a.sqrtrem())
(xmpz(11), xmpz(2))
>>> print(b.sqrtrem())
(xmpz(21), xmpz(15))
>>> for i in range(5):
...    a=_g.xmpz(123)
...    b=_g.xmpz(456)
...    print(a.root(i+1),b.root(i+1))
...
(xmpz(123), 1) (xmpz(456), 1)
(xmpz(11), 0) (xmpz(21), 0)
(xmpz(4), 0) (xmpz(7), 0)
(xmpz(3), 0) (xmpz(4), 0)
(xmpz(2), 0) (xmpz(3), 0)
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
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
xmpz(3)
>>> temp=_g.gcdext(a,b)
>>> temp[0]==a*temp[1]+b*temp[2]
True
>>> _g.lcm(a,b)
xmpz(18696)
>>> _g.set_prefer_mutable(1)
>>> _g.fac(7)
xmpz(5040)
>>> _g.fib(17)
xmpz(1597)
>>> _g.set_prefer_mutable(0)
>>> del temp
>>> for i in range(10):
...     print(_g.bincoef(10,i))
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
>>> _g.divm(b,a,_g.xmpz(20))
xmpz(12)
>>> _g.divm(a,b,100)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.divm(a,b,100)
ZeroDivisionError: not invertible
>>> _g.set_prefer_mutable(1)
>>> _g.divm(6,12,14)
xmpz(4)
>>> _g.divm(0,1,2)
xmpz(0)
>>> _g.set_prefer_mutable(0)
>>> a.invert(100)
>>> a
xmpz(87)
>>> b.invert(100)
>>> b
xmpz(0)
>>> _g.invert(3)
Traceback (innermost last):
  ...
TypeError: invert() requires 'mpz','mpz' arguments
>>> a.invert()
Traceback (innermost last):
  ...
TypeError: invert() takes exactly one argument (0 given)
>>> a.invert(3, 4)
Traceback (innermost last):
  ...
TypeError: invert() takes exactly one argument (2 given)
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
ValueError: factor must be > 0
>>> _g.remove(3)
Traceback (innermost last):
  ...
TypeError: remove() requires 'mpz','mpz' arguments
>>> a.remove()
Traceback (innermost last):
  ...
TypeError: remove() requires 'mpz','mpz' arguments
>>> a.remove(3, 4)
Traceback (innermost last):
  ...
TypeError: remove() requires 'mpz','mpz' arguments
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
TypeError: jacobi() requires 'mpz','mpz' arguments
>>> a.jacobi()
Traceback (innermost last):
  ...
TypeError: jacobi() requires 'mpz','mpz' arguments
>>> a.jacobi(3, 4)
Traceback (innermost last):
  ...
TypeError: jacobi() requires 'mpz','mpz' arguments
>>> _g.legendre(10,3)
1
>>> _g.legendre(10,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: legendre's y must be odd and > 0
>>> _g.legendre(3)
Traceback (innermost last):
  ...
TypeError: legendre() requires 'mpz','mpz' arguments
>>> a.legendre()
Traceback (innermost last):
  ...
TypeError: legendre() requires 'mpz','mpz' arguments
>>> a.legendre(3, 4)
Traceback (innermost last):
  ...
TypeError: legendre() requires 'mpz','mpz' arguments
>>> _g.kronecker(10,3)
1
>>> _g.kronecker(10,-3)
1
>>> _g.kronecker(3)
Traceback (innermost last):
  ...
TypeError: kronecker() requires 'mpz','mpz' arguments
>>> a.kronecker()
Traceback (innermost last):
  ...
TypeError: kronecker() requires 'mpz','mpz' arguments
>>> a.kronecker(3, 4)
Traceback (innermost last):
  ...
TypeError: kronecker() requires 'mpz','mpz' arguments
>>> a=10**20
>>> b=a+39
>>> _g.jacobi(a,b)
1
>>> _g.legendre(a,b)
1
>>> _g.kronecker(a,b)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: Either arg in Kronecker must fit in an int
'''

def _test(chat=None):
    if chat:
        print("Unit tests for gmpy2 (mpz functionality)")
        print("    running on Python %s" % sys.version)
        print()
        if _g.gmp_version():
            print("Testing gmpy2 %s (GMP %s), default caching (%s, %s)" % (
                (_g.version(), _g.gmp_version(), _g.get_cache()[0],
                _g.get_cache()[1])))
        else:
            print("Testing gmpy2 %s (MPIR %s), default caching (%s, %s)" % (
                (_g.version(), _g.mpir_version(), _g.get_cache()[0],
                _g.get_cache()[1])))
    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat: print("Repeating tests, with caching disabled")
    _g.set_cache(0,128)

    sav = sys.stdout
    class _Dummy:
        encoding = None
        def write(self,*whatever):
            pass
    try:
        sys.stdout = _Dummy()
        doctest.testmod(thismod, report=0)
    finally:
        sys.stdout = sav

    if chat:
        print()
        print("Overall results for mpz:")
    return doctest.master.summarize(chat)


if __name__=='__main__':
    _test(1)


