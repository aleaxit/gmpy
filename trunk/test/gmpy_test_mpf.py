# partial unit test for gmpy2 mpf functionality
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> filter(lambda x: not x.startswith('__'), dir(a))
['_copy', 'acos', 'acosh', 'add', 'agm', 'ai', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 'binary', 'cbrt', 'ceil', 'cos', 'cosh', 'cot', 'coth', 'csc', 'csch', 'digamma', 'digits', 'div', 'eint', 'erf', 'erfc', 'exp', 'exp10', 'exp2', 'expm1', 'f2q', 'floor', 'gamma', 'hypot', 'is_inf', 'is_nan', 'is_number', 'is_regular', 'is_zero', 'j0', 'j1', 'jn', 'lgamma', 'li2', 'lngamma', 'log', 'log10', 'log1p', 'log2', 'max', 'min', 'mul', 'next_above', 'next_below', 'next_toward', 'pow', 'precision', 'qdiv', 'rec_sqrt', 'reldiff', 'root', 'round', 'sec', 'sech', 'sign', 'sin', 'sin_cos', 'sinh', 'sinh_cosh', 'sqrt', 'square', 'sub', 'tan', 'tanh', 'trunc', 'y0', 'y1', 'yn', 'zeta']
>>>
'''
import sys

import gmpy2 as _g, doctest, sys
__test__={}
a=_g.mpf('123.456')
b=_g.mpf('789.123')

__test__['functions']=\
r'''
>>> _g.log(2)
mpf('6.9314718055994529e-1')
>>> _g.log(10)
mpf('2.3025850929940459e0')
>>> _g.log('a')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: log() requires 'mpf' argument
>>> _g.log(float('nan'))
mpf('nan')
>>> _g.log(float('inf'))
mpf('inf')
>>> _g.log(float('-inf'))
mpf('nan')
>>> _g.mpf('12.3456').log()
mpf('2.5132997242892183e0')
>>> _g.mpf('12.3456').log(7)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: log() takes no arguments (1 given)
'''

__test__['elemop']=\
r'''
>>> str(a+b)
'912.57900000000006'
>>> str(a-b)
'-665.66700000000003'
>>> str(a*b)
'97421.969088000013'
>>> str(a/b)
'0.15644709379906555'
>>> str(b+a)
'912.57900000000006'
>>> str(b-a)
'665.66700000000003'
>>> str(b*a)
'97421.969088000013'
>>> str(b/a)
'6.3919372083981338'
>>> str(-a)
'-123.456'
>>> _g.mpf(2,200) + 3
mpf('5.0e0')
>>> 3 + _g.mpf(2,200)
mpf('5.0e0')
>>> _g.mpf(2,200) * 3
mpf('6.0e0')
>>> 3 * _g.mpf(2,200)
mpf('6.0e0')
>>> str(abs(-a))
'123.456'
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
>>> import math
>>> math.ceil(a)
124.0
>>> str(a.ceil())
'124.0'
>>> str(_g.ceil(a))
'124.0'
>>> math.floor(a)
123.0
>>> str(a.floor())
'123.0'
>>> str(_g.floor(a))
'123.0'
>>> str(a.trunc())
'123.0'
>>> str(_g.trunc(a))
'123.0'
>>> x=-a
>>> math.floor(x)
-124.0
>>> str(x.floor())
'-124.0'
>>> str(_g.floor(x))
'-124.0'
>>> str(x.ceil())
'-123.0'
>>> math.ceil(x)
-123.0
>>> str(_g.ceil(x))
'-123.0'
>>> str(x.trunc())
'-123.0'
>>> str(_g.trunc(x))
'-123.0'
>>> _g.ceil(12.3)==math.ceil(12.3)
1
>>> _g.floor(12.3)==math.floor(12.3)
1
>>> _g.ceil(-12.3)==math.ceil(-12.3)
1
>>> _g.floor(-12.3)==math.floor(-12.3)
1
>>> (a**2).reldiff(float(a)**2) < 1.03 * (2.0**-(a.precision-1))
1
>>> (a**2).reldiff(a*a) < (2.0**-(a.precision-1))
1
>>> (b**2).reldiff(float(b)**2) < 1.03 * (2.0**-(b.precision-1))
1
>>> (b**2).reldiff(b*b) < (2.0**-(b.precision-1))
1
>>> _g.reldiff(3.4)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff(3.4)
TypeError: reldiff() requires 'mpf,'mpf' arguments
>>> a.reldiff()
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff()
TypeError: reldiff() requires 'mpf,'mpf' arguments
>>> a.reldiff(3, 4)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff(3, 4)
TypeError: reldiff() requires 'mpf,'mpf' arguments
>>> a.sqrt()
mpf('1.1111075555498667e1')
>>> _g.sqrt(a)
mpf('1.1111075555498667e1')
>>> _g.sqrt(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
  File "a.py", line 9, in _er
    raise ValueError, what
ValueError: sqrt() of negative number
>>> _g.pi()
mpf('3.1415926535897931e0')
>>> _g.pi(64)
mpf('3.14159265358979323851e0',64)
>>> _g.pi(200)
mpf('3.1415926535897932384626433832795028841971693993751058209749445e0',200)
>>> save=_g.get_precision()
>>> _g.const_catalan()
mpf('9.1596559417721901e-1')
>>> _g.const_euler()
mpf('5.7721566490153287e-1')
>>> _g.const_log2()
mpf('6.9314718055994529e-1')
>>> _g.const_pi()
mpf('3.1415926535897931e0')
>>> _g.set_precision(100)
>>> _g.const_catalan()
mpf('9.1596559417721901505460351493252e-1',100)
>>> _g.const_euler()
mpf('5.7721566490153286060651209008234e-1',100)
>>> _g.const_log2()
mpf('6.9314718055994530941723212145798e-1',100)
>>> _g.const_pi()
mpf('3.1415926535897932384626433832793e0',100)
>>> _g.set_precision(save)
>>> del(save)
>>> import pickle
>>> flt = _g.mpf(1234.6789)
>>> flt == pickle.loads(pickle.dumps(flt))
True
>>> flt = _g.mpf('1.1')
>>> flt == pickle.loads(pickle.dumps(flt))
True
'''


from gmpy_truediv import truediv
__test__['newdiv']=\
r'''
>>>
>>> a/b
mpf('1.5644709379906555e-1')
>>> a//b
mpf('0.0e0')
>>> truediv(a,b)
mpf('1.5644709379906555e-1')
>>> b/a
mpf('6.3919372083981338e0')
>>> b//a
mpf('6.0e0')
>>> truediv(b,a)
mpf('6.3919372083981338e0')
>>>
'''

__test__['cmpr']=\
r'''
>>> c=_g.mpf(a)
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
>>> not _g.mpf(0)
1
>>> not a
0
>>> a.f2q(0.1)
mpz(123)
>>> a.f2q(0.01)
mpz(123)
>>> a.f2q(0.001)
mpq(247,2)
>>> a.f2q(0.0001)
mpq(1358,11)
>>> a.f2q(0.00001)
mpq(7037,57)
>>> a.f2q(0.000001)
mpq(15432,125)
>>> a.f2q(0.0000001)
mpq(15432,125)
>>> a.f2q()
mpq(15432,125)
>>> print _g.mpf(_g.mpz(1234))
1234.0
>>> x=1000*1000*1000*1000L
>>> _g.mpf(x)
mpf('1.0e12')
>>> c=_g.mpf(a)
>>> a is c
1
>>> c=_g.mpf(a,99)
>>> a is c
0
>>> a==c
1
>>> _g.mpf('1.1') == _g.mpf('1.1') * _g.mpf(1)
True
>>> _g.mpf('1.1',64) == _g.mpf('1.1',128)
False
>>> _g.mpf('1.1',64) == _g.mpf(_g.mpf('1.1',128),64)
True
>>> a = _g.mpf('.123', 64)
>>> b = _g.mpf('.123', 128)
>>> c = _g.mpf('.123', 128) * _g.mpf(1, 128)
>>> a == b
False
>>> a == c
False
>>> b == c
False
>>> a == b.round(64)
True
>>> a == _g.round(b, 64)
True
>>> _g.mpf('ffffffffffffffffe8000000000000000', 256, 16).round(64).digits(16)
('10000000000000000', 34, 64)
>>> _g.mpf('fffffffffffffffff8000000000000000', 256, 16).round(64).digits(16)
('10000000000000000', 34, 64)
>>> b.round(64)
mpf('1.22999999999999999999e-1',64)
'''

__test__['format']=\
r'''
>>> str(a)
'123.456'
>>> repr(a)
"mpf('1.23456e2')"
>>> a.digits(10,0)
('12345600000000000', 3, 53)
>>> a.digits(10,1)
Traceback (most recent call last):
  ...
ValueError: digits must be 0 or >= 2
>>> a.digits(10,2)
('12', 3, 53)
>>> a.digits(10,3)
('123', 3, 53)
>>> a.digits(10,4)
('1235', 3, 53)
>>> a.digits(10,5)
('12346', 3, 53)
>>> a.digits(10,6)
('123456', 3, 53)
>>> a.digits(10,7)
('1234560', 3, 53)
>>> a.digits(10,8)
('12345600', 3, 53)
>>> junk=_g.set_fcoform(14)
>>> frmt=_g.set_fcoform(14)
>>> frmt
'%.14e'
>>> ofmt=_g.set_fcoform(frmt)
>>> ofmt
'%.14e'
>>> _g.mpf(3.4)
mpf('3.3999999999999999e0')
>>> print _g.mpf(3.4)
3.3999999999999999
>>> _g.set_fcoform(junk)
'%.14e'
>>> a.digits(1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be in the interval 2 ... 62
>>> a.digits(2,-1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: digits must be 0 or >= 2
>>> a.digits(10,0,0,-1,2)
Traceback (most recent call last):
  ...
TypeError: function takes at most 2 arguments (5 given)
>>> saveprec=a.precision
>>> newa = a.round(33)
>>> newa
mpf('1.23456e2',33)
>>> newa = newa.round(saveprec)
>>> newa.precision==saveprec
1
>>> del(newa)
>>> _g.digits(_g.mpf(23.45))
('23449999999999999', 2, 53)
>>> _g.binary('pep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: binary() requires a gmpy2 object as argument
>>>
'''

__test__['binio']=\
r'''
>>> epsilon=_g.mpf(2)**-(a.precision)
>>> ba=a.binary()
>>> a.reldiff(_g.mpf(ba,0,256)) <= epsilon
1
>>> len(ba)
16
>>> for i in range(len(ba)):
...     print ord(ba[i]),
...     if i==len(ba)-1: print
...
8 53 0 0 0 1 0 0 0 123 116 188 106 126 249 220
>>> na=(-a).binary()
>>> (-a).reldiff(_g.mpf(na,0,256)) <= epsilon
1
>>> na[0] == chr(ord(ba[0])|1)
1
>>> for bd,nd in zip(ba[1:],na[1:]):
...    assert bd==nd
>>> ia=(1/a).binary()
>>> (1/a).reldiff(_g.mpf(ia,0,256)) <= epsilon
1
>>> _g.binary(_g.mpf(0))
'\x04'
>>> _g.mpf(_g.binary(_g.mpf(0)), 0, 256) == 0
1
>>> _g.binary(_g.mpf(0.5))
'\x085\x00\x00\x00\x00\x00\x00\x00\x80'
>>> _g.mpf(_g.binary(_g.mpf(0.5)), 0, 256) == 0.5
1
>>> _g.binary(_g.mpf(-0.5))
'\t5\x00\x00\x00\x00\x00\x00\x00\x80'
>>> _g.mpf(_g.binary(_g.mpf(-0.5)), 0, 256) == -0.5
1
>>> _g.binary(_g.mpf(-2.0))
'\t5\x00\x00\x00\x01\x00\x00\x00\x02'
>>> _g.mpf(_g.binary(_g.mpf(-2.0)), 0, 256) == -2.0
1
>>> _g.binary(_g.mpf(2.0))
'\x085\x00\x00\x00\x01\x00\x00\x00\x02'
>>> _g.mpf(_g.binary(_g.mpf(2.0)), 0, 256) == 2.0
1
>>> prec=_g.set_precision(0)
Traceback (most recent call last):
  ...
ValueError: invalid value for precision
>>> junk=_g.set_precision(53)
>>> hash(_g.mpf(23.0))==hash(23)
1
>>> print _g.mpf('\004',0,256)
0.0e0
>>> long(a)
123L
>>> long(-a)
-123L
>>> int(a)
123
>>> int(-a)
-123
>>>
'''

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy2 (mpf functionality)"
        print "    running on Python %s" % sys.version
        print
        if _g.gmp_version():
            print "Testing gmpy2 %s (GMP %s) with default caching (%s, %s)" % (
                (_g.version(), _g.gmp_version(), _g.get_cache()[0],
                _g.get_cache()[1]))
        else:
            print "Testing gmpy2 %s (MPIR %s) with default caching (%s, %s)" % (
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
        print "Overall results for mpf:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)

