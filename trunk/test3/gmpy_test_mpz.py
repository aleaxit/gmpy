# partial unit test for gmpy2 mpz functionality
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> list(filter(lambda x: not x.startswith('_'), dir(_g)))
['ModeMPFR', 'ModePython', 'RoundAwayZero', 'RoundDown', 'RoundToNearest', 'RoundToZero', 'RoundUp', 'acos', 'acosh', 'add', 'agm', 'ai', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 'binary', 'bincoef', 'bit_clear', 'bit_flip', 'bit_length', 'bit_mask', 'bit_scan0', 'bit_scan1', 'bit_set', 'bit_test', 'cbrt', 'cdiv', 'cdiv2exp', 'cdivmod', 'cdivmod2exp', 'ceil', 'clear_erangeflag', 'clear_flags', 'clear_inexactflag', 'clear_nanflag', 'clear_overflow', 'clear_underflow', 'cmod', 'cmod2exp', 'comb', 'const_catalan', 'const_euler', 'const_log2', 'const_pi', 'cos', 'cosh', 'cot', 'coth', 'csc', 'csch', 'denom', 'digamma', 'digits', 'div', 'divexact', 'divm', 'eint', 'erf', 'erfc', 'exp', 'exp10', 'exp2', 'expm1', 'f2q', 'fac', 'factorial', 'fdiv', 'fdiv2exp', 'fdivmod', 'fdivmod2exp', 'fib', 'fib2', 'floor', 'fma', 'fmod', 'fmod2exp', 'fms', 'gamma', 'gcd', 'gcdext', 'get_cache', 'get_emax', 'get_emax_max', 'get_emin', 'get_emin_min', 'get_max_precision', 'get_mode', 'get_precision', 'get_rounding', 'get_ternary', 'hamdist', 'hypot', 'inf', 'invert', 'is_erangeflag', 'is_even', 'is_inexactflag', 'is_inf', 'is_lessgreater', 'is_nan', 'is_nanflag', 'is_number', 'is_odd', 'is_overflow', 'is_power', 'is_prime', 'is_regular', 'is_square', 'is_underflow', 'is_unordered', 'is_zero', 'j0', 'j1', 'jacobi', 'jn', 'kronecker', 'lcm', 'legendre', 'lgamma', 'li2', 'license', 'lngamma', 'log', 'log10', 'log1p', 'log2', 'lucas', 'lucas2', 'max', 'min', 'mp_limbsize', 'mp_version', 'mpf', 'mpfr_version', 'mpq', 'mpz', 'mul', 'nan', 'next_above', 'next_below', 'next_prime', 'next_toward', 'numdigits', 'numer', 'pack', 'pi', 'popcount', 'pow', 'qdiv', 'rec_sqrt', 'reldiff', 'remove', 'root', 'rootrem', 'round', 'sec', 'sech', 'set_cache', 'set_debug', 'set_emax', 'set_emin', 'set_erangeflag', 'set_inexactflag', 'set_mode', 'set_nanflag', 'set_overflow', 'set_precision', 'set_rounding', 'set_underflow', 'sign', 'sin', 'sin_cos', 'sinh', 'sinh_cosh', 'sqrt', 'sqrtrem', 'square', 'sub', 'tan', 'tanh', 'tdiv', 'tdiv2exp', 'tdivmod', 'tdivmod2exp', 'tmod', 'tmod2exp', 'trunc', 'unpack', 'version', 'xbit_mask', 'xmpz', 'y0', 'y1', 'yn', 'zero', 'zeta']
>>> dir(a)
['__abs__', '__add__', '__and__', '__bool__', '__class__', '__delattr__', '__divmod__', '__doc__', '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__iadd__', '__ifloordiv__', '__ilshift__', '__imod__', '__imul__', '__index__', '__init__', '__int__', '__invert__', '__ipow__', '__irshift__', '__isub__', '__le__', '__len__', '__lshift__', '__lt__', '__mod__', '__mul__', '__ne__', '__neg__', '__new__', '__or__', '__pos__', '__pow__', '__radd__', '__rand__', '__rdivmod__', '__reduce__', '__reduce_ex__', '__repr__', '__rfloordiv__', '__rlshift__', '__rmod__', '__rmul__', '__ror__', '__rpow__', '__rrshift__', '__rshift__', '__rsub__', '__rtruediv__', '__rxor__', '__setattr__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__xor__', '_copy', 'binary', 'bincoef', 'bit_clear', 'bit_flip', 'bit_length', 'bit_scan0', 'bit_scan1', 'bit_set', 'bit_test', 'cdiv', 'cdiv2exp', 'cmod', 'cmod2exp', 'comb', 'digits', 'divexact', 'fdiv', 'fdiv2exp', 'fmod', 'fmod2exp', 'hamdist', 'invert', 'is_even', 'is_odd', 'is_power', 'is_prime', 'is_square', 'jacobi', 'kronecker', 'legendre', 'next_prime', 'numdigits', 'popcount', 'qdiv', 'remove', 'root', 'rootrem', 'sign', 'sqrt', 'sqrtrem', 'tdiv', 'tdiv2exp', 'tmod', 'tmod2exp']
>>>
'''
import gmpy2 as _g, doctest, sys, operator, gc
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
    from functools import reduce
    assert reduce(mul, factors)==savex
    return factors

__test__['index']=\
r'''
>>> range(333)[a]
123
>>> range(333)[b]
Traceback (innermost last):
  ...
IndexError: range object index out of range
'''

__test__['elemop']=\
r'''
>>> a+b
mpz(579)
>>> a-b
mpz(-333)
>>> a*b
mpz(56088)
>>> a//b
mpz(0)
>>> a/b
mpf('2.6973684210526316e-1')
>>> b//a
mpz(3)
>>> b/a
mpf('3.7073170731707319e0')
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
>>> (-a).sign()
-1
>>> z=b-b; z.sign()
0
>>> import pickle
>>> pickle.loads(pickle.dumps(_g.mpz(12346789)))
mpz(12346789)
>>> s='12345678901234567890123456789'
>>> int(s) == _g.mpz(s)
True
>>> _g.mpz(s) == int(s)
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
>>> _g.get_mode() == _g.ModePython
True
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
>>> -a * float('Inf')
mpf('-inf')
>>> float('Inf') * -a
mpf('-inf')
>>> -a * float('-Inf')
mpf('inf')
>>> float('-Inf') * -a
mpf('inf')
>>> a * float('nan')
mpf('nan')
>>> float('nan') * a
mpf('nan')
>>> _g.mpz(0) * float('Inf')
mpf('nan')
>>> _g.mpz(0) * float('-Inf')
mpf('nan')
>>> float('Inf') * _g.mpz(0)
mpf('nan')
>>> float('-Inf') * _g.mpz(0)
mpf('nan')
>>> a / float('Inf')
mpf('0.0e0')
>>> -a / float('Inf')
mpf('-0.0e0')
>>> float('Inf') / a
mpf('inf')
>>> float('Inf') / -a
mpf('-inf')
>>> a / float('-Inf')
mpf('-0.0e0')
>>> -a / float('-Inf')
mpf('0.0e0')
>>> float('-Inf') / a
mpf('-inf')
>>> float('-Inf') / -a
mpf('inf')
>>> a / float('nan')
mpf('nan')
>>> float('nan') / a
mpf('nan')
>>> float('nan') / _g.mpz(0)
Traceback (most recent call last):
  ...
ZeroDivisionError: mpf division by zero
>>> divmod(a, float('Inf'))
(mpf('0.0e0'), mpf('1.23e2'))
>>> divmod(a, float('-Inf'))
(mpf('-1.0e0'), mpf('-inf'))
>>> divmod(-a, float('Inf'))
(mpf('-1.0e0'), mpf('inf'))
>>> divmod(-a, float('-Inf'))
(mpf('0.0e0'), mpf('-1.23e2'))
>>> divmod(a, float('nan'))
(mpf('nan'), mpf('nan'))
>>> divmod(-a, float('nan'))
(mpf('nan'), mpf('nan'))
>>> divmod(_g.mpz(0), float('Inf'))
(mpf('0.0e0'), mpf('0.0e0'))
>>> divmod(_g.mpz(0), float('-Inf'))
(mpf('-0.0e0'), mpf('-0.0e0'))
>>> divmod(_g.mpz(0), float('nan'))
(mpf('nan'), mpf('nan'))
>>> divmod(float('Inf'), a)
(mpf('nan'), mpf('nan'))
>>> divmod(float('-Inf'), a)
(mpf('nan'), mpf('nan'))
>>> divmod(float('Inf'), -a)
(mpf('nan'), mpf('nan'))
>>> divmod(float('-Inf'), -a)
(mpf('nan'), mpf('nan'))
>>> divmod(float('nan'), a)
(mpf('nan'), mpf('nan'))
>>> divmod(float('nan'), -a)
(mpf('nan'), mpf('nan'))
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

__test__['divexact']=\
r'''
>>> a=_g.mpz('1234567912345678912345679')
>>> b=_g.mpz('789789789789789789789789')
>>> c=a*b
>>> _g.divexact(c,a)
mpz(789789789789789789789789)
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
>>> _g.mpz(23)
mpz(23)
>>> _g.mpz(-23)
mpz(-23)
>>> x=1000*1000*1000*1000*1000*1000*1000
>>> _g.mpz(x)
mpz(1000000000000000000000)
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
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.hamdist()
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.hamdist(3, 4)
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> a.fmod2exp(5)
mpz(27)
>>> b.fmod2exp(5)
mpz(8)
>>> b.fmod2exp(5)==(b%32)
1
>>> a.fmod2exp(5)==(a%32)
1
>>> a.bit_set(20)
mpz(1048699)
>>> a.bit_clear(0)
mpz(122)
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
>>> n=_g.mpz(-(7+6*16+5*256+7*4092))
>>> [n.bit_scan0(j) for j in range(18)]
[1, 1, 3, 3, 6, 6, 6, 8, 8, 10, 10, 12, 12, 13, 14, -1, None, None]
>>> [n.bit_scan1(j) for j in range(33)]
[0, 2, 2, 4, 4, 5, 7, 7, 9, 9, 11, 11, 15, 15, 15, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
>>> _g.mpz(0).bit_length()
0
>>> _g.mpz(12345).bit_length()
14
>>> _g.bit_mask(9)
mpz(511)
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
'0o173'
>>> _g.mpz('123')
mpz(123)
>>> _g.mpz('1001001011',2)
mpz(587)
>>> _g.mpz('1001001011',2).digits(2)
'0b1001001011'
>>> for i in range(2,63):
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
3C
39
36
33
30
2d
2b
2Z
2X
2V
2T
2R
2P
2N
2L
2J
2H
2F
2D
2B
29
27
25
23
21
1z
>>> print(a.digits(63))
Traceback (innermost last):
  ...
ValueError: base must be either 0 or in the interval 2 ... 62
>>> _g.mpz('43')
mpz(43)
>>> _g.mpz('043')
mpz(43)
>>> _g.mpz('43',0)
mpz(43)
>>> _g.mpz('0o43')
mpz(35)
>>> _g.mpz('0x43')
mpz(67)
>>> _g.mpz('0x43',10)
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
b'{'
>>> _g.mpz(ba,256)
mpz(123)
>>> _g.mpz(ba,256)==a
1
>>> _g.binary(_g.mpz(123))
b'{'
>>> z=_g.mpz('melancholy',256)
>>> z
mpz(573406620562849222387053)
>>> int(z)
573406620562849222387053
>>> divmod(z,a)
(mpz(4661842443600400182008), mpz(69))
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
ValueError: base must be either 0 or in the interval 2 ... 62
>>> hash(a)
123
>>> hash(b)
456
>>> hash(z) == hash(int(z))
True
>>> _g.mpz(_g.binary(_g.mpz(-123)),256)
mpz(-123)
>>> int(_g.mpz(-3))
-3
'''

__test__['number']=\
r'''
>>> a.rootrem(2)
(mpz(11), mpz(2))
>>> a.rootrem(3)
(mpz(4), mpz(59))
>>> _g.rootrem(a*a)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: rootrem() requires 'mpz','int' arguments
>>> _g.rootrem(a*a,2)
(mpz(123), mpz(0))
>>> print(a.sqrt())
11
>>> print(b.sqrt())
21
>>> print(a.sqrtrem())
(mpz(11), mpz(2))
>>> print(b.sqrtrem())
(mpz(21), mpz(15))
>>> for i in range(5):
...    print(a.root(i+1),b.root(i+1))
...
(mpz(123), True) (mpz(456), True)
(mpz(11), False) (mpz(21), False)
(mpz(4), False) (mpz(7), False)
(mpz(3), False) (mpz(4), False)
(mpz(2), False) (mpz(3), False)
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
>>> for x in range(45678):
...     _xx=_g.divm(b,a,20)
>>> del _xx
>>> __ = gc.collect()
>>> (_memsize()-_siz) <= 32
True
>>> a.invert(100)
mpz(87)
>>> b.invert(100)
mpz(0)
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
ValueError: sqrt() of negative number
>>> _g.sqrtrem(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: sqrtrem() of negative number
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
        print("    on Python %s" % sys.version)
        print("Testing gmpy2 {}".format(_g.version()))
        print("  Mutliple-precision library:   {}".format(_g.mp_version()))
        print("  Floating-point library:       {}".format(_g.mpfr_version()))
        print("  Caching Values: (Number)      {}".format(_g.get_cache()[0]))
        print("  Caching Values: (Size, limbs) {}".format(_g.get_cache()[1]))
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


