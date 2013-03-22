# partial unit test for gmpy2 mpz functionality
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> filter(lambda x: not x.startswith('_'), dir(_g))
['Default', 'DivisionByZeroError', 'ExponentOutOfBoundsError', 'InexactResultError', 'InvalidOperationError', 'OverflowResultError', 'RangeError', 'RoundAwayZero', 'RoundDown', 'RoundToNearest', 'RoundToZero', 'RoundUp', 'UnderflowResultError', 'acos', 'acosh', 'add', 'agm', 'ai', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 'bincoef', 'bit_clear', 'bit_flip', 'bit_length', 'bit_mask', 'bit_scan0', 'bit_scan1', 'bit_set', 'bit_test', 'c_div', 'c_div_2exp', 'c_divmod', 'c_divmod_2exp', 'c_mod', 'c_mod_2exp', 'cbrt', 'ceil', 'check_range', 'comb', 'const_catalan', 'const_euler', 'const_log2', 'const_pi', 'context', 'copy_sign', 'cos', 'cosh', 'cot', 'coth', 'csc', 'csch', 'degrees', 'denom', 'digamma', 'digits', 'div', 'div_2exp', 'divexact', 'divm', 'eint', 'erf', 'erfc', 'exp', 'exp10', 'exp2', 'expm1', 'f2q', 'f_div', 'f_div_2exp', 'f_divmod', 'f_divmod_2exp', 'f_mod', 'f_mod_2exp', 'fac', 'factorial', 'fib', 'fib2', 'floor', 'fma', 'fmod', 'fms', 'frac', 'frexp', 'from_binary', 'fsum', 'gamma', 'gcd', 'gcdext', 'get_cache', 'get_context', 'get_emax_max', 'get_emin_min', 'get_exp', 'get_max_precision', 'hamdist', 'hypot', 'ieee', 'inf', 'invert', 'iroot', 'iroot_rem', 'is_bpsw_prp', 'is_euler_prp', 'is_even', 'is_extra_strong_lucas_prp', 'is_fermat_prp', 'is_fibonacci_prp', 'is_finite', 'is_inf', 'is_infinite', 'is_integer', 'is_lessgreater', 'is_lucas_prp', 'is_nan', 'is_number', 'is_odd', 'is_power', 'is_prime', 'is_regular', 'is_selfridge_prp', 'is_signed', 'is_square', 'is_strong_bpsw_prp', 'is_strong_lucas_prp', 'is_strong_prp', 'is_strong_selfridge_prp', 'is_unordered', 'is_zero', 'isqrt', 'isqrt_rem', 'j0', 'j1', 'jacobi', 'jn', 'kronecker', 'lcm', 'legendre', 'lgamma', 'li2', 'license', 'lngamma', 'local_context', 'log', 'log10', 'log1p', 'log2', 'lucas', 'lucas2', 'lucasu', 'lucasu_mod', 'lucasv', 'lucasv_mod', 'max2', 'maxnum', 'min2', 'minnum', 'modf', 'mp_limbsize', 'mp_version', 'mpc', 'mpc_random', 'mpc_version', 'mpfr', 'mpfr_from_old_binary', 'mpfr_grandom', 'mpfr_random', 'mpfr_version', 'mpq', 'mpq_from_old_binary', 'mpz', 'mpz_from_old_binary', 'mpz_random', 'mpz_rrandomb', 'mpz_urandomb', 'mul', 'mul_2exp', 'nan', 'next_above', 'next_below', 'next_prime', 'next_toward', 'norm', 'num_digits', 'numer', 'pack', 'phase', 'polar', 'popcount', 'powmod', 'proj', 'qdiv', 'radians', 'random_state', 'rec_sqrt', 'rect', 'reldiff', 'remainder', 'remove', 'remquo', 'rint', 'rint_ceil', 'rint_floor', 'rint_round', 'rint_trunc', 'root', 'round2', 'round_away', 'sec', 'sech', 'set_cache', 'set_context', 'set_exp', 'set_sign', 'sign', 'sin', 'sin_cos', 'sinh', 'sinh_cosh', 'sqrt', 'square', 'sub', 't_div', 't_div_2exp', 't_divmod', 't_divmod_2exp', 't_mod', 't_mod_2exp', 'tan', 'tanh', 'to_binary', 'trunc', 'unpack', 'version', 'xbit_mask', 'xmpz', 'y0', 'y1', 'yn', 'zero', 'zeta']
>>> filter(lambda x: not x.startswith('__'), dir(a))
['bit_clear', 'bit_flip', 'bit_length', 'bit_scan0', 'bit_scan1', 'bit_set', 'bit_test', 'denominator', 'digits', 'num_digits', 'numerator']
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
    import gmpy2 as _g
    savex=x
    prime=2
    x=_g.mpz(x)
    factors=[]
    while x>=prime:
        newx,mult=_g.remove(x,prime)
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
>>> _g.sign(a)
1
>>> _g.sign(-a)
-1
>>> z=b-b; _g.sign(z)
0
>>> _g.mpz(4)**(0.5)
mpfr('2.0')
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
>>> _g.is_even(2)
True
'''

from gmpy_truediv import truediv
__test__['newdiv']=\
r'''
>>> a/b
mpz(0)
>>> a//b
mpz(0)
>>> truediv(a,b)
mpfr('0.26973684210526316')
>>> b/a
mpz(3)
>>> b//a
mpz(3)
>>> truediv(b,a)
mpfr('3.7073170731707319')
>>>
'''

__test__['divexact']=\
r'''
>>> a=_g.mpz('1234567912345678912345679')
>>> b=_g.mpz('789789789789789789789789')
>>> c=a*b
>>> _g.divexact(c,a)
mpz(789789789789789789789789L)
>>>
'''

__test__['divmod']=\
r'''
>>> _g.c_divmod(17,5)
(mpz(4), mpz(-3))
>>> _g.c_divmod(-17,5)
(mpz(-3), mpz(-2))
>>> _g.c_divmod(17,-5)
(mpz(-3), mpz(2))
>>> _g.c_divmod(-17,-5)
(mpz(4), mpz(3))
>>> _g.f_divmod(17,5)
(mpz(3), mpz(2))
>>> _g.f_divmod(-17,5)
(mpz(-4), mpz(3))
>>> _g.f_divmod(17,-5)
(mpz(-4), mpz(-3))
>>> _g.f_divmod(-17,-5)
(mpz(3), mpz(-2))
>>> _g.t_divmod(17,5)
(mpz(3), mpz(2))
>>> _g.t_divmod(-17,5)
(mpz(-3), mpz(-2))
>>> _g.t_divmod(17,-5)
(mpz(-3), mpz(2))
>>> _g.t_divmod(-17,-5)
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
>>> _g.context().trap_divzero == False
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
mpfr('inf')
>>> float('Inf') + a
mpfr('inf')
>>> a + float('-Inf')
mpfr('-inf')
>>> float('-Inf') + a
mpfr('-inf')
>>> a + float('nan')
mpfr('nan')
>>> float('nan') + a
mpfr('nan')
>>> a - float('Inf')
mpfr('-inf')
>>> float('Inf') - a
mpfr('inf')
>>> a - float('-Inf')
mpfr('inf')
>>> float('-Inf') - a
mpfr('-inf')
>>> a - float('nan')
mpfr('nan')
>>> float('nan') - a
mpfr('nan')
>>> a * float('Inf')
mpfr('inf')
>>> float('Inf') * a
mpfr('inf')
>>> a * float('-Inf')
mpfr('-inf')
>>> float('-Inf') * a
mpfr('-inf')
>>> -a * float('Inf')
mpfr('-inf')
>>> float('Inf') * -a
mpfr('-inf')
>>> -a * float('-Inf')
mpfr('inf')
>>> float('-Inf') * -a
mpfr('inf')
>>> a * float('nan')
mpfr('nan')
>>> float('nan') * a
mpfr('nan')
>>> _g.mpz(0) * float('Inf')
mpfr('nan')
>>> _g.mpz(0) * float('-Inf')
mpfr('nan')
>>> float('Inf') * _g.mpz(0)
mpfr('nan')
>>> float('-Inf') * _g.mpz(0)
mpfr('nan')
>>> a / float('Inf')
mpfr('0.0')
>>> -a / float('Inf')
mpfr('-0.0')
>>> float('Inf') / a
mpfr('inf')
>>> float('Inf') / -a
mpfr('-inf')
>>> a / float('-Inf')
mpfr('-0.0')
>>> -a / float('-Inf')
mpfr('0.0')
>>> float('-Inf') / a
mpfr('-inf')
>>> float('-Inf') / -a
mpfr('inf')
>>> a / float('nan')
mpfr('nan')
>>> float('nan') / a
mpfr('nan')
>>> float('nan') / _g.mpz(0)
mpfr('nan')
>>> float('nan') / _g.mpz(0)
mpfr('nan')
>>> divmod(a, float('Inf'))
(mpfr('0.0'), mpfr('123.0'))
>>> divmod(a, float('-Inf'))
(mpfr('-1.0'), mpfr('-inf'))
>>> divmod(-a, float('Inf'))
(mpfr('-1.0'), mpfr('inf'))
>>> divmod(-a, float('-Inf'))
(mpfr('0.0'), mpfr('-123.0'))
>>> divmod(a, float('nan'))
(mpfr('nan'), mpfr('nan'))
>>> divmod(-a, float('nan'))
(mpfr('nan'), mpfr('nan'))
>>> divmod(_g.mpz(0), float('Inf'))
(mpfr('0.0'), mpfr('0.0'))
>>> divmod(_g.mpz(0), float('-Inf'))
(mpfr('-0.0'), mpfr('-0.0'))
>>> divmod(_g.mpz(0), float('nan'))
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('Inf'), a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('-Inf'), a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('Inf'), -a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('-Inf'), -a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('nan'), a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('nan'), -a)
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('Inf'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('-Inf'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
>>> divmod(float('nan'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
>>> divmod(_g.mpfr('Inf'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
>>> divmod(_g.mpfr('-Inf'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
>>> divmod(_g.mpfr('nan'), _g.mpz(0))
(mpfr('nan'), mpfr('nan'))
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
>>> _g.popcount(a)
6
>>> _g.popcount(b)
4
>>> _g.popcount(-7)
-1
>>> _g.popcount(0)
0
>>> _g.hamdist(a,b)
6
>>> _g.hamdist(3)
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> _g.hamdist()
Traceback (innermost last):
  ...
TypeError: hamdist() requires 'mpz','mpz' arguments
>>> _g.f_mod_2exp(a,5)
mpz(27)
>>> _g.f_mod_2exp(b,5)
mpz(8)
>>> _g.f_mod_2exp(b,5)==(b%32)
1
>>> _g.f_mod_2exp(a,5)==(a%32)
1
>>> a.bit_set(20)
mpz(1048699)
>>> a.bit_clear(0)
mpz(122)
>>> for i in range(8):
...     print a.bit_test(i),
...     if i==7: print
...
True True False True True True True False
>>> for i in range(10):
...     print b.bit_test(i),
...     if i==9: print
...
False False False True False False True True True False
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
  ...
ValueError: base must be in the interval 2 ... 62
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
>>> a== _g.from_binary(_g.to_binary(a))
True
>>> -a== _g.from_binary(_g.to_binary(-a))
True
>>> z=_g.mpz(573406620562849222387053L)
>>> divmod(z,a)
(mpz(4661842443600400182008L), mpz(69))
>>> for i in range(2,37):
...    print i,z.num_digits(i),
...    if i%6==0: print
...
2 79 3 50 4 40 5 35 6 31
7 29 8 27 9 25 10 24 11 23 12 23
13 22 14 21 15 21 16 20 17 20 18 19
19 19 20 19 21 18 22 18 23 18 24 18
25 18 26 17 27 17 28 17 29 17 30 17
31 16 32 16 33 16 34 16 35 16 36 16
>>> _g.num_digits(23)
2
>>> _g.num_digits(23,2)
5
>>> _g.num_digits(23,99)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be in the interval 2 ... 62
>>> hash(a)
123
>>> hash(b)
456
>>> hash(z) == hash(long(z))
True
>>> long(_g.mpz(-3))
-3L
>>>
'''

__test__['number']=\
r'''
>>> _g.iroot_rem(a,2)
(mpz(11), mpz(2))
>>> _g.iroot_rem(a,3)
(mpz(4), mpz(59))
>>> _g.iroot_rem(a*a)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: iroot_rem() requires 'mpz','int' arguments
>>> _g.iroot_rem(a*a,2)
(mpz(123), mpz(0))
>>> print _g.isqrt(a), _g.isqrt(b)
11 21
>>> print _g.isqrt_rem(a), _g.isqrt_rem(b)
(mpz(11), mpz(2)) (mpz(21), mpz(15))
>>> for i in range(5):
...    print _g.iroot(a,i+1),_g.iroot(b,i+1)
...
(mpz(123), True) (mpz(456), True)
(mpz(11), False) (mpz(21), False)
(mpz(4), False) (mpz(7), False)
(mpz(3), False) (mpz(4), False)
(mpz(2), False) (mpz(3), False)
>>> _g.is_square(a)
0
>>> _g.is_power(a)
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
>>> _g.invert(a,100)
mpz(87)
>>> _g.invert(b,100)
mpz(0)
>>> _g.invert(3)
Traceback (innermost last):
  ...
TypeError: invert() requires 'mpz','mpz' arguments
>>> _g.comb(3,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: binomial coefficient with negative k
>>> _g.isqrt(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: isqrt() of negative number
>>> _g.isqrt_rem(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: isqrt_rem() of negative number
>>> _g.remove(3,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: factor must be > 1
>>> _g.remove(3)
Traceback (innermost last):
  ...
TypeError: remove() requires 'mpz','mpz' arguments
>>> _g.is_prime(3,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: repetition count for is_prime() must be positive
>>> _g.jacobi(10,3)
1
>>> _g.jacobi(10,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: y must be odd and >0
>>> _g.jacobi(3)
Traceback (innermost last):
  ...
TypeError: jacobi() requires 'mpz','mpz' arguments
>>> _g.legendre(10,3)
1
>>> _g.legendre(10,-3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: y must be odd and >0
>>> _g.legendre(3)
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
>>> a=10L**20
>>> b=a+39
>>> _g.jacobi(a,b)
1
>>> _g.legendre(a,b)
1
>>> _g.kronecker(a,b)
1
'''

def _test(chat=None):
    try:
        float('nan')
    except ValueError:
        if "special" in __test__:
            del(__test__["special"])

    if chat:
        print "Unit tests for gmpy2 (mpz functionality)"
        print "    on Python %s" % sys.version
        print "Testing gmpy2 {0}".format(_g.version())
        print "  Mutliple-precision library:   {0}".format(_g.mp_version())
        print "  Floating-point library:       {0}".format(_g.mpfr_version())
        print "  Complex library:              {0}".format(_g.mpc_version())
        print "  Caching Values: (Number)      {0}".format(_g.get_cache()[0])
        print "  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1])

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

