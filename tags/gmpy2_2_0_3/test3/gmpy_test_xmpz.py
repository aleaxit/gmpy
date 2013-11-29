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
mpfr('0.26973684210526316')
>>> b//a
mpz(3)
>>> b/a
mpfr('3.7073170731707319')
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
>>> _g.sign(a)
1
>>> -a
>>> _g.sign(a)
-1
>>> a=_g.xmpz(123)
>>> z=b-b; _g.sign(z)
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
>>> -a
>>> a * float('Inf')
mpfr('-inf')
>>> float('Inf') * a
mpfr('-inf')
>>> a * float('-Inf')
mpfr('inf')
>>> float('-Inf') * a
mpfr('inf')
>>> -a
>>> a * float('nan')
mpfr('nan')
>>> float('nan') * a
mpfr('nan')
>>> _g.xmpz(0) * float('Inf')
mpfr('nan')
>>> _g.xmpz(0) * float('-Inf')
mpfr('nan')
>>> float('Inf') * _g.xmpz(0)
mpfr('nan')
>>> float('-Inf') * _g.xmpz(0)
mpfr('nan')
>>> a / float('Inf')
mpfr('0.0')
>>> float('Inf') / a
mpfr('inf')
>>> a / float('-Inf')
mpfr('-0.0')
>>> float('-Inf') / a
mpfr('-inf')
>>> a / float('nan')
mpfr('nan')
>>> float('nan') / a
mpfr('nan')
>>> -a
>>> a / float('Inf')
mpfr('-0.0')
>>> float('Inf') / a
mpfr('-inf')
>>> a / float('-Inf')
mpfr('0.0')
>>> float('-Inf') / a
mpfr('inf')
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
>>> d=a.copy()
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
>>> a.bit_set(20)
mpz(1048699)
>>> a=_g.xmpz(123)
>>> a.bit_clear(0)
mpz(122)
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
>>> _g.bit_mask(9)
mpz(511)
>>> _g.xbit_mask(9)
xmpz(511)
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
ValueError: base must be in the interval 2 ... 62
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

__test__['number']=\
r'''
>>> a=_g.xmpz(123)
>>> b=_g.xmpz(456)
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
'''

def _test(chat=None):
    if chat:
        print("Unit tests for gmpy2 (xmpz functionality)")
        print("    on Python %s" % sys.version)
        print("Testing gmpy2 {0}".format(_g.version()))
        print("  Mutliple-precision library:   {0}".format(_g.mp_version()))
        print("  Floating-point library:       {0}".format(_g.mpfr_version()))
        print("  Complex library:              {0}".format(_g.mpc_version()))
        print("  Caching Values: (Number)      {0}".format(_g.get_cache()[0]))
        print("  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1]))
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


