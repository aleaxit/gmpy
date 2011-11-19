# partial unit test for gmpy mpf functionality
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> filter(lambda x: not x.startswith('__'), dir(a))
['_copy', 'binary', 'ceil', 'digits', 'f2q', 'floor', 'getprec', 'getrprec', 'qdiv', 'reldiff', 'round', 'setprec', 'sign', 'sqrt', 'trunc']
>>>
'''
import warnings
warnings.filterwarnings('ignore', 'setprec')

import sys

import gmpy as _g, doctest, sys
__test__={}
a=_g.mpf('123.456')
b=_g.mpf('789.123')

__test__['elemop']=\
r'''
>>> str(a+b)
'912.579'
>>> str(a-b)
'-665.667'
>>> str(a*b)
'97421.969088'
>>> str(a/b)
'0.156447093799065544915'
>>> str(b+a)
'912.579'
>>> str(b-a)
'665.667'
>>> str(b*a)
'97421.969088'
>>> str(b/a)
'6.39193720839813374807'
>>> str(-a)
'-123.456'
>>> _g.mpf(2,200) + 3
mpf('5.e0',200)
>>> 3 + _g.mpf(2,200)
mpf('5.e0',200)
>>> _g.mpf(2,200) * 3
mpf('6.e0',200)
>>> 3 * _g.mpf(2,200)
mpf('6.e0',200)
>>> str(abs(-a))
'123.456'
>>> _g.fsign(b-a)
1
>>> _g.fsign(b-b)
0
>>> _g.fsign(a-b)
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
>>> (a**2).reldiff(float(a)**2) < 1.03 * (2.0**-(a.getrprec()-1))
1
>>> (a**2).reldiff(a*a) < (2.0**-(a.getprec()-1))
1
>>> (b**2).reldiff(float(b)**2) < 1.03 * (2.0**-(b.getrprec()-1))
1
>>> (b**2).reldiff(b*b) < (2.0**-(b.getprec()-1))
1
>>> _g.reldiff(3.4)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff(3.4)
TypeError: function takes exactly 2 arguments (1 given)
>>> a.reldiff()
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff()
TypeError: function takes exactly 1 argument (0 given)
>>> a.reldiff(3, 4)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff(3, 4)
TypeError: function takes exactly 1 argument (2 given)
>>> a.sqrt()
mpf('1.11110755554986664846e1')
>>> _g.fsqrt(a)
mpf('1.11110755554986664846e1')
>>> _g.fsqrt(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
  File "a.py", line 9, in _er
    raise ValueError, what
ValueError: sqrt of negative number
>>> _g.pi(64)
mpf('3.14159265358979323846e0',64)
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
mpf('1.56447093799065544915e-1')
>>> a//b
mpf('0.e0')
>>> truediv(a,b)
mpf('1.56447093799065544915e-1')
>>> b/a
mpf('6.39193720839813374807e0')
>>> b//a
mpf('6.e0')
>>> truediv(b,a)
mpf('6.39193720839813374807e0')
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
mpf('1.e12')
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
True
>>> a == b.round(64)
True
>>> a == _g.fround(b, 64)
True
>>> _g.mpf('ffffffffffffffffe8000000000000000', 256, 16).round(64).digits(16)
'f.fffffffffffffffe@32'
>>> _g.mpf('fffffffffffffffff8000000000000000', 256, 16).round(64).digits(16)
'1.@33'
>>> b.round(64)
mpf('1.23e-1',64)
'''

__test__['format']=\
r'''
>>> str(a)
'123.456'
>>> repr(a)
"mpf('1.23456e2')"
>>> _g.set_tagoff(0)
1
>>> a
gmpy.mpf('1.23456e2')
>>> _g.set_tagoff(1)
0
>>> a.digits(10,0)
'1.23456e2'
>>> a.digits(10,1)
'1.e2'
>>> a.digits(10,2)
'1.2e2'
>>> a.digits(10,3)
'1.23e2'
>>> a.digits(10,4)
'1.235e2'
>>> a.digits(10,5)
'1.2346e2'
>>> a.digits(10,6)
'1.23456e2'
>>> a.digits(10,7)
'1.23456e2'
>>> a.digits(10,8)
'1.23456e2'
>>> for i in range(11,99):
...     tempa='%.16f' % (i/10.0)
...     tempb=_g.mpf(i/10.0).digits(10,17)
...     assert tempb.startswith(tempa.rstrip('0')), (tempa, tempb)
...
>>> junk=_g.set_fcoform(14)
>>> frmt=_g.set_fcoform(14)
>>> frmt
'%.14e'
>>> ofmt=_g.set_fcoform(frmt)
>>> ofmt
'%.14e'
>>> _g.mpf(3.4)
mpf('3.39999999999999999998e0')
>>> print _g.mpf(3.4)
3.39999999999999999998
>>> _g.set_fcoform(junk)
'%.14e'
>>> a.digits(1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be either 0 or in the interval 2 ... 62
>>> a.digits(2,-1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: digits must be >= 0
>>> a.digits(10,0,0,-1,2)
('123456', 3, 53)
>>> saveprec=a.getrprec()
>>> a.setprec(33)
>>> a
mpf('1.23456e2',33)
>>> a.setprec(saveprec)
>>> a.getrprec()==saveprec
1
>>> _g.fdigits(2.2e5, 0, 6, -10, 10)
'220000.0'
>>> _g.fdigits(2.2e-5, 0, 6, -10, 10)
'0.000022'
>>> _g.digits(_g.mpf(23.45))
Traceback (most recent call last):
  ...
TypeError: digits() expects 'mpz',['int'] arguments
>>> _g.fbinary('pep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: argument can not be converted to mpf
>>>
'''

__test__['binio']=\
r'''
>>> epsilon=_g.mpf(2)**-(a.getrprec())
>>> ba=a.binary()
>>> a.reldiff(_g.mpf(ba,0,256)) <= epsilon
1
>>> len(ba)
18
>>> for i in range(len(ba)):
...     print ord(ba[i]),
...     if i==len(ba)-1: print
...
8 53 0 0 0 1 0 0 0 123 116 188 106 126 249 219 34 209
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
>>> _g.fbinary(0)
'\x04'
>>> _g.mpf(_g.fbinary(0), 0, 256) == 0
1
>>> _g.fbinary(0.5)
'\x085\x00\x00\x00\x00\x00\x00\x00\x80'
>>> _g.mpf(_g.fbinary(0.5), 0, 256) == 0.5
1
>>> _g.fbinary(-0.5)
'\t5\x00\x00\x00\x00\x00\x00\x00\x80'
>>> _g.mpf(_g.fbinary(-0.5), 0, 256) == -0.5
1
>>> _g.fbinary(-2.0)
'\t5\x00\x00\x00\x01\x00\x00\x00\x02'
>>> _g.mpf(_g.fbinary(-2.0), 0, 256) == -2.0
1
>>> _g.fbinary(2.0)
'\x085\x00\x00\x00\x01\x00\x00\x00\x02'
>>> _g.mpf(_g.fbinary(2.0), 0, 256) == 2.0
1
>>> prec=_g.set_minprec(0)
>>> junk=_g.set_minprec(prec)
>>> a.getrprec()==prec
1
>>> b.getrprec()==prec
1
>>> _g.mpf(1.0).getrprec()==prec
1
>>> hash(_g.mpf(23.0))==hash(23)
1
>>> print _g.mpf('\004',0,256)
0.0
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
        print "Unit tests for gmpy 1.15 (mpf functionality)"
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
        print "Overall results for mpf:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)

