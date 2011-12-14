# partial unit test for gmpy extra cover
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> print int(_g.gmp_version()[:3] in ('5.0', '4.3', ''))
1
>>> print int(_g.mpir_version()[:3] in ('2.3', '2.4', '2.5', ''))
1
>>> _g.version()
'1.15'
>>> int('gmpy.c' in _g._cvsid())
1
'''

import gmpy as _g, doctest, sys
__test__={}
r = _g.rand

__test__['misc_stuff']=\
r'''
>>> junk=_g.set_debug(0)
>>> knuj=_g.set_debug(junk)
>>> junk==_g.set_debug(junk)
1
>>> _g.set_fcoform(None)
>>> _g.set_fcoform()
>>> print _g.mpf(3.0)
3.0
>>> _g.gmp_limbsize() in (32, 64)
True
>>> _g.mpz(u"123")
mpz(123)
>>> _g.mpq(u"12/37")
mpq(12,37)
>>> _g.mpf(u"123")
mpf('1.23e2')
'''

try:
    x = float('inf')
except ValueError:
    pass
else:
    __test__['infinity'] = \
r'''
>>> x = float('inf')
>>> n = float('nan')
>>> _g.mpf(x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpf(-x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpq(x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpq(-x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpz(x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpz(-x)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle infinity
>>> _g.mpf(n)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle nan
>>> _g.mpf(n)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle nan
>>> _g.mpf(n)
Traceback (most recent call last):
  ...
ValueError: gmpy does not handle nan
'''

__test__['user_errors']=\
r'''
>>> _g.version(23)
Traceback (most recent call last):
  ...
TypeError: version expects 0 arguments
>>> _g.gmp_version(23)
Traceback (most recent call last):
  ...
TypeError: gmp_version expects 0 arguments
>>> _g.get_cache(23)
Traceback (most recent call last):
  ...
TypeError: get_cache expects 0 arguments
>>> _g.set_cache()
Traceback (most recent call last):
  ...
TypeError: function takes exactly 2 arguments (0 given)
>>> _g.set_cache(200)
Traceback (most recent call last):
  ...
TypeError: function takes exactly 2 arguments (1 given)
>>> _g.set_cache(200,-23)
Traceback (most recent call last):
  ...
ValueError: object size must between 0 and 16384
>>> _g.set_cache(2000,256)
Traceback (most recent call last):
  ...
ValueError: cache must between 0 and 1000
>>> _g.set_cache(-23,256)
Traceback (most recent call last):
  ...
ValueError: cache must between 0 and 1000
>>> _g.set_cache(200,256000)
Traceback (most recent call last):
  ...
ValueError: object size must between 0 and 16384
>>> _g.set_debug()
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes exactly 1 argument (0 given)
>>> _g.set_debug(2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes exactly 1 argument (2 given)
>>> _g.set_debug('boh')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.set_minprec(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: minimum precision must be >= 0
>>> _g.set_fcoform(33)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: number of digits n must be 0<n<=30
>>> _g.set_fcoform([])
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: set_fcoform argument must be int, string, or None
>>> _g.mpz('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string without NULL characters expected
>>> _g.mpf('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string without NULL characters expected
>>> _g.mpq('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string without NULL characters expected
>>> _g.mpq('bo',256)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (too short)
>>> _g.mpq('bologna',256)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (num len)
>>> _g.mpq('\001\000\000\000\003\002',256)
mpq(3,2)
>>> _g.mpq('\002\000\000\000\003\377\002',256)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (num sgn)
>>> _g.mpq('\001\000\000\000\003\002\377',256)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (den sgn)
>>> _g.mpq('ba/bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid digits
>>> print 'ba/bo'
ba/bo
>>> _g.mpq('1/bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid digits
>>> print '1/bo'
1/bo
>>> _g.mpq('1/0')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: mpq: zero denominator
>>> _g.mpf([])
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpf() expects numeric or string argument
>>> _g.mpf('bo',0,256)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string too short to be a gmpy.mpf binary encoding
>>> _g.mpf('bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid digits
>>> int(_g.mpz(1000L*1000*1000*1000*1000*1000*1000))
1000000000000000000000L
>>> _g.scan0(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: starting bit must be >= 0
>>> _g.scan1(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: starting bit must be >= 0
>>> _g.lowbits(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: nbits must be > 0
>>> _g.getbit(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: bit_index must be >= 0
>>> _g.setbit(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: bit_index must be >= 0
>>> _g.mpz(23).setbit(12,1,2,3)
Traceback (most recent call last):
  ...
TypeError: setbit() expects 'mpz','int'[,'int'] arguments
>>> _g.setbit(12,1,2,3)
Traceback (most recent call last):
  ...
TypeError: setbit() expects 'mpz','int'[,'int'] arguments
>>> _g.root(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: n must be > 0
>>> _g.root(12,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: n must be > 0
>>> _g.root(-12,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: root of negative number
>>> _g.digits(3.14)
Traceback (most recent call last):
  ...
TypeError: digits() expects 'mpz',['int'] arguments
>>> _g.digits(3,'peep')
Traceback (most recent call last):
  ...
TypeError: digits() expects 'mpz',['int'] arguments
>>> _g.fdigits(3.14,'peep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.qdigits(3.14,'peep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: argument can not be converted to mpq
>>> _g.qdigits(3,'peep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.mpz(3).digits('bu')
Traceback (most recent call last):
  ...
TypeError: digits() expects 'mpz',['int'] arguments
>>> _g.mpf(3).digits('bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.qdiv(3,_g.mpq(1))
mpz(3)
>>> _g.qdiv(3,_g.mpz(1))
mpz(3)
>>> _g.qdiv(3,_g.mpf(1))
mpz(3)
>>> _g.qdiv(3,1.0)
mpz(3)
>>> _g.qdiv(3,1L)
mpz(3)
>>> _g.qdiv(3)
mpz(3)
>>> _g.qdiv(3,'bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: second argument can not be converted to mpq
>>> _g.mpq(2).qdiv(3,4)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes at most 1 argument (2 given)
>>> _g.qdiv(3,4,5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes at most 2 arguments (3 given)
>>> _g.qdiv('bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: first argument can not be converted to mpq
>>> _g.qdiv(1.0,1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: first argument can not be converted to mpq
>>> _g.qdiv(1,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: qdiv: zero divisor
>>> _g.f2q(-1.0)
mpz(-1)
>>> _g.mpz(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpz() requires 1 or 2 arguments
>>> _g.mpz('bi','bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpz(): base must be an integer
>>> _g.mpz('bi',99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for gmpy.mpz must be 0, 256, or in the interval 2 ... 62 .
>>> _g.mpz(1,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpz() with numeric argument needs exactly 1 argument
>>> _g.mpz(None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpz() expects numeric or string argument
>>> _g.mpq(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpq() requires 1 or 2 arguments
>>> _g.mpq('bi','bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpq(): base must be an integer
>>> _g.mpq('bi',99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for gmpy.mpq() must be 0, 256, or in the interval 2 ... 36 .
>>> _g.mpq(None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpq() expects numeric or string argument
>>> _g.mpq(1,None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: argument can not be converted to mpq
>>> _g.mpq(1,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: mpq: zero denominator
>>> _g.mpf()
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpf() requires 1 to 3 arguments
>>> _g.mpf(1,'bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpf(): bits must be an integer
>>> _g.mpf(1,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: bits for gmpy.mpf must be >= 0
>>> _g.mpf('ba',0,'bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpf(): base must be an integer
>>> _g.mpf('ba',0,99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for gmpy.mpf must be 0, 256, or in the interval 2 ... 62 .
>>> _g.mpf(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: gmpy.mpf() with numeric 1st argument needs 1 or 2 arguments
>>> +_g.mpz(1)
mpz(1)
>>> +_g.mpf(1)
mpf('1.e0')
>>> +_g.mpq(1)
mpq(1,1)
>>> _g.mpz(2)**-2
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpz.pow with negative power
>>> _g.mpz(2)**_g.mpz(1000000*10000000000000L)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpz.pow outrageous exponent
>>> pow(_g.mpz(2),3,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpz.pow divide by zero
>>> pow(_g.mpz(2),3,-5)
mpz(-2)
>>> pow(_g.mpq(2),3,-5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow no modulo allowed
>>> a=10000000000L**2
>>> _g.mpq(2)**a
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow outrageous exp num
>>> _g.mpq(2)**_g.mpq(1,a)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow outrageous exp den
>>> _g.mpq(2)**0
mpq(1,1)
>>> _g.mpq(2)**-1
mpq(1,2)
>>> _g.mpq(2)**_g.mpq(1,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow fractional exponent, inexact-root
>>> _g.mpq(-2)**_g.mpq(1,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow fractional exponent, nonreal-root
>>> _g.mpq(0)**_g.mpq(1,2)
mpq(0,1)
>>> _g.mpq(0)**-1
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: mpq.pow 0 base to <0 exponent
>>> _g.mpq(-1)**-1
mpq(-1,1)
>>> _g.mpf(9,100)**2
mpf('8.1e1',100)
>>> _g.mpf(9,100)**0.5
mpf('3.e0',100)
>>> _g.mpf(9,100)**_g.mpf(0.5)
mpf('3.e0')
>>> _g.mpf(0)**2
mpf('0.e0')
>>> pow(_g.mpf(2),3,-5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpf.pow no modulo allowed
>>> _g.mpz(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpz' and 'str'
>>> _g.mpz(1)+_g.mpf(1)
mpf('2.e0')
>>> _g.mpz(1)+_g.mpq(1)
mpq(2,1)
>>> _g.mpq(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpq' and 'str'
>>> _g.mpf(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpf' and 'str'
>>> _g.mpf(1)+_g.mpq(2)
mpf('3.e0')
>>> divmod(_g.mpz(3),0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: mpz divmod by zero
>>> divmod(_g.mpz(0),3)
(mpz(0), mpz(0))
>>> _g.divm(1,2,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: not invertible
>>> abs(_g.mpq(0))
mpq(0,1)
>>> _g.mpz(0)**2
mpz(0)
>>> _g.mpq(-2)**0
mpq(1,1)
>>> _g.fac(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: factorial of negative number
>>> _g.fib(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: Fibonacci of negative number
>>> _g.comb(7,-3)
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
>>> _g.fsqrt(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: sqrt of negative number
>>> _g.jacobi(23, -34)
Traceback (most recent call last):
  ...
ValueError: jacobi's y must be odd prime > 0
>>> _g.legendre(23, -34)
Traceback (most recent call last):
  ...
ValueError: legendre's y must be odd and > 0
>>> # guard against conversion error on 64-bit systems
>>> _g.mpz(2**32) != _g.mpz(0)
True
>>> # test hash properties on 64-bit systems
>>> temp = 123456789012345678901234567890
>>> hash(temp) == hash(_g.mpz(temp))
True
>>> del temp
'''

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy 1.15 (extra cover)"
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

    if chat:
        print
        print "Overall results for cvr:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)
