# partial unit test for gmpy2 mpfr functionality
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> list([x for x in dir(a) if x != '__dir__'])
['__abs__', '__add__', '__bool__', '__ceil__', '__class__', '__delattr__', '__divmod__', '__doc__', '__eq__', '__float__', '__floor__', '__floordiv__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__int__', '__le__', '__lt__', '__mod__', '__mul__', '__ne__', '__neg__', '__new__', '__pos__', '__pow__', '__radd__', '__rdivmod__', '__reduce__', '__reduce_ex__', '__repr__', '__rfloordiv__', '__rmod__', '__rmul__', '__round__', '__rpow__', '__rsub__', '__rtruediv__', '__setattr__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__trunc__', 'as_integer_ratio', 'as_mantissa_exp', 'as_simple_fraction', 'conjugate', 'digits', 'imag', 'is_integer', 'precision', 'rc', 'real']
>>>
'''
import sys

import gmpy2 as _g, doctest, sys
__test__={}
a=_g.mpfr('123.456')
b=_g.mpfr('789.123')

__test__['functions']=\
r'''
>>> _g.log(2)
mpfr('0.69314718055994529')
>>> _g.log(10)
mpfr('2.3025850929940459')
>>> _g.log('a')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: log() argument type not supported
>>> _g.log(float('nan'))
mpfr('nan')
>>> _g.log(float('inf'))
mpfr('inf')
>>> _g.log(float('-inf'))
mpfr('nan')
>>> _g.log(_g.mpfr('12.3456'))
mpfr('2.5132997242892183')
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
>>> str(abs(-a))
'123.456'
>>> _g.mpfr(2) + 3
mpfr('5.0')
>>> 3 + _g.mpfr(2)
mpfr('5.0')
>>> _g.mpfr(2) * 3
mpfr('6.0')
>>> 3 * _g.mpfr(2)
mpfr('6.0')
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
>>> import math
>>> math.ceil(a)
mpfr('124.0')
>>> str(a.__ceil__())
'124.0'
>>> str(_g.ceil(a))
'124.0'
>>> math.floor(a)
mpfr('123.0')
>>> str(a.__floor__())
'123.0'
>>> str(_g.floor(a))
'123.0'
>>> math.trunc(a)
mpfr('123.0')
>>> str(a.__trunc__())
'123.0'
>>> str(_g.trunc(a))
'123.0'
>>> x=-a
>>> math.floor(x)
mpfr('-124.0')
>>> str(x.__floor__())
'-124.0'
>>> str(_g.floor(x))
'-124.0'
>>> math.ceil(x)
mpfr('-123.0')
>>> str(x.__ceil__())
'-123.0'
>>> str(_g.ceil(x))
'-123.0'
>>> math.trunc(x)
mpfr('-123.0')
>>> str(x.__trunc__())
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
>>> _g.reldiff(a**2,float(a)**2) < 1.03 * (2.0**-(a.precision-1))
1
>>> _g.reldiff(a**2,a*a) < (2.0**-(a.precision-1))
1
>>> _g.reldiff(b**2,float(b)**2) < 1.03 * (2.0**-(b.precision-1))
1
>>> _g.reldiff(b**2,b*b) < (2.0**-(b.precision-1))
1
>>> _g.reldiff(3.4)
Traceback (innermost last):
  File "<pyshell#184>", line 1, in ?
    _g.reldiff(3.4)
TypeError: reldiff() requires 'mpfr','mpfr' arguments
>>> _g.sqrt(a)
mpfr('11.111075555498667')
>>> _g.sqrt(-1)
mpfr('nan')
>>> save=_g.get_context().precision
>>> _g.const_catalan()
mpfr('0.91596559417721901')
>>> _g.const_euler()
mpfr('0.57721566490153287')
>>> _g.const_log2()
mpfr('0.69314718055994529')
>>> _g.const_pi()
mpfr('3.1415926535897931')
>>> _g.get_context().precision = 100
>>> _g.const_catalan()
mpfr('0.91596559417721901505460351493252',100)
>>> _g.const_euler()
mpfr('0.57721566490153286060651209008234',100)
>>> _g.const_log2()
mpfr('0.69314718055994530941723212145798',100)
>>> _g.const_pi()
mpfr('3.1415926535897932384626433832793',100)
>>> _g.get_context().precision = save
>>> del(save)
>>> import pickle
>>> flt = _g.mpfr(1234.6789)
>>> flt == pickle.loads(pickle.dumps(flt))
True
>>> flt = _g.mpfr('1.1')
>>> flt == pickle.loads(pickle.dumps(flt))
True
'''

__test__['newdiv']=\
r'''
>>>
>>> a/b
mpfr('0.15644709379906555')
>>> a//b
mpfr('0.0')
>>> b/a
mpfr('6.3919372083981338')
>>> b//a
mpfr('6.0')
>>>
'''

__test__['cmpr']=\
r'''
>>> c=_g.mpfr(a)
>>> c is a
1
>>> c==a
1
>>> c>a
0
>>> c<a
0
>>> a>b
0
>>> a<b
1
>>> not _g.mpfr(0)
1
>>> not a
0
>>> _g.f2q(a,0.1)
mpz(123)
>>> _g.f2q(a,0.01)
mpz(123)
>>> _g.f2q(a,0.001)
mpq(247,2)
>>> _g.f2q(a,0.0001)
mpq(1358,11)
>>> _g.f2q(a,0.00001)
mpq(7037,57)
>>> _g.f2q(a,0.000001)
mpq(15432,125)
>>> _g.f2q(a,0.0000001)
mpq(15432,125)
>>> _g.f2q(a)
mpq(15432,125)
>>> print(_g.mpfr(_g.mpz(1234)))
1234.0
>>> x=1000*1000*1000*1000
>>> _g.mpfr(x)
mpfr('1000000000000.0')
>>> c=_g.mpfr(a)
>>> a is c
1
>>> c=_g.mpfr(a,99)
>>> a is c
0
>>> a==c
1
>>> _g.mpfr('1.1') == _g.mpfr('1.1') * _g.mpfr(1)
True
>>> _g.mpfr('1.1',64) == _g.mpfr('1.1',128)
False
>>> _g.mpfr('1.1',64) == _g.mpfr(_g.mpfr('1.1',128),64)
True
>>> a = _g.mpfr('.123', 64)
>>> b = _g.mpfr('.123', 128)
>>> c = _g.mpfr('.123', 128) * _g.mpfr(1, 128)
>>> a == b
False
>>> a == c
False
>>> b == c
False
>>> f=_g.mpfr
>>> bool(float('nan')) == bool(f('nan'))
True
>>> (float('nan') == float('nan')) == (float('nan') == f('nan')) == (f('nan') == float('nan')) == (f('nan') == f('nan'))
True
>>> (float('nan') != float('nan')) == (float('nan') != f('nan')) == (f('nan') != float('nan')) == (f('nan') != f('nan'))
True
>>> bool(float('inf')) == bool(f('inf'))
True
>>> bool(float('-inf')) == bool(f('-inf'))
True
'''

__test__['format']=\
r'''
>>> str(a)
'123.456'
>>> repr(a)
"mpfr('123.456')"
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
>>> for i in range(11,99):
...     tempa=('%.16f' % (i/10.0)).replace('.','')
...     tempb=_g.mpfr(i/10.0).digits(10,17)[0]
...     assert tempb.startswith(tempa.rstrip('0')), (tempa, tempb)
...
>>> _g.mpfr(3.4)
mpfr('3.3999999999999999')
>>> print(_g.mpfr(3.4))
3.3999999999999999
>>> a.digits(1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: base must be in the interval 2 ... 62
>>> a.digits(2,-1)
Traceback (most recent call last):
  File "<string>", line 1, in ?
ValueError: digits must be 0 or >= 2
>>> a.digits(10,0)
('12345600000000000', 3, 53)
>>> saveprec=a.precision
>>> newa = _g.round2(a,33)
>>> newa
mpfr('123.456',33)
>>> newa = _g.round2(newa,saveprec)
>>> newa.precision==saveprec
1
>>> del(newa)
>>> _g.digits(_g.mpfr(23.45))
('23449999999999999', 2, 53)
>>>
'''

__test__['binio']=\
r'''
>>> _g.from_binary('pep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: from_binary() requires bytes argument
>>> a == _g.from_binary(_g.to_binary(a))
True
>>> -a == _g.from_binary(_g.to_binary(-a))
True
>>> 1/a == _g.from_binary(_g.to_binary(1/a))
True
>>> _g.to_binary(_g.mpfr(0))
b'\x04\x00\x00\x005\x00\x00\x00'
>>> 0 == _g.from_binary(_g.to_binary(_g.mpfr(0)))
True
>>> 0.5 == _g.from_binary(_g.to_binary(_g.mpfr(0.5)))
True
>>> -0.5 == _g.from_binary(_g.to_binary(_g.mpfr(-0.5)))
True
>>> -2.0 == _g.from_binary(_g.to_binary(_g.mpfr(-2.0)))
True
>>> 2.0 == _g.from_binary(_g.to_binary(_g.mpfr(2.0)))
True
>>> hash(_g.mpfr(23.0))==hash(23)
True
>>> int(a)
123
>>> int(-a)
-123
>>>
'''

def _test(chat=None):
    if chat:
        print("Unit tests for gmpy2 (mpfr functionality)")
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
        doctest.testmod(thismod, report=0, optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
    finally:
        sys.stdout = sav

    if chat:
        print()
        print("Overall results for mpfr:")
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)

