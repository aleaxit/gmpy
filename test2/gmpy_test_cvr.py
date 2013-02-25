# partial unit test for gmpy2 extra cover
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> _g.version()
'2.0.0b5'
>>> int('gmpy2.c' in _g._cvsid())
1
'''

import gmpy2 as _g, doctest, sys
__test__={}

__test__['misc_stuff']=\
r'''
>>> print _g.mpfr(3.0)
3.0
>>> _g.mp_limbsize() in (32, 64)
True
>>> _g.mpz(u"123")
mpz(123)
>>> _g.mpq(u"12/37")
mpq(12,37)
>>> _g.mpfr(u"123")
mpfr('123.0')
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
>>> _g.mpfr(x)
mpfr('inf')
>>> _g.mpfr(-x)
mpfr('-inf')
>>> _g.mpq(x)
Traceback (most recent call last):
  ...
OverflowError: 'mpq' does not support Infinity
>>> _g.mpq(-x)
Traceback (most recent call last):
  ...
OverflowError: 'mpq' does not support Infinity
>>> _g.mpz(x)
Traceback (most recent call last):
  ...
OverflowError: 'mpz' does not support Infinity
>>> _g.mpz(-x)
Traceback (most recent call last):
  ...
OverflowError: 'mpz' does not support Infinity
>>> _g.mpfr(n)
mpfr('nan')
>>> _g.mpq(n)
Traceback (most recent call last):
  ...
ValueError: 'mpq' does not support NaN
>>> _g.mpz(n)
Traceback (most recent call last):
  ...
ValueError: 'mpz' does not support NaN
'''

__test__['user_errors']=\
r'''
>>> _g.version(23)
Traceback (most recent call last):
  ...
TypeError: version() takes no arguments (1 given)
>>> _g.mp_version(23)
Traceback (most recent call last):
  ...
TypeError: mp_version() takes no arguments (1 given)
>>> _g.get_cache(23)
Traceback (most recent call last):
  ...
TypeError: get_cache() takes no arguments (1 given)
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
ValueError: cache size must between 0 and 1000
>>> _g.set_cache(-23,256)
Traceback (most recent call last):
  ...
ValueError: cache size must between 0 and 1000
>>> _g.set_cache(200,256000)
Traceback (most recent call last):
  ...
ValueError: object size must between 0 and 16384
>>> _g.context().precision = -1
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid value for precision
>>> _g.mpz('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string contains NULL characters
>>> _g.mpfr('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid digits
>>> _g.mpq('12'+chr(0)+'34')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: string contains NULL characters
>>> _g.mpq_from_old_binary('bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (too short)
>>> _g.mpq_from_old_binary('bologna')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpq binary (num len)
>>> _g.mpq_from_old_binary('\001\000\000\000\003\002')
mpq(3,2)
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
ZeroDivisionError: zero denominator in 'mpq'
>>> _g.mpfr([])
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpfr() requires numeric or string argument
>>> _g.mpfr_from_old_binary('bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid mpf binary encoding (too short)
>>> _g.mpfr('bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: invalid digits
>>> int(_g.mpz(1000L*1000*1000*1000*1000*1000*1000))
1000000000000000000000L
>>> _g.bit_scan0(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: starting bit must be >= 0
>>> _g.bit_scan1(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: starting bit must be >= 0
>>> _g.bit_test(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: bit_index must be >= 0
>>> _g.bit_set(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: bit_index must be >= 0
>>> _g.mpz(23).bit_set(12,1,2,3)
Traceback (most recent call last):
  ...
TypeError: bit_set() takes exactly one argument (4 given)
>>> _g.bit_set(12,1,2,3)
Traceback (most recent call last):
  ...
TypeError: bit_set() requires 'mpz','int' arguments
>>> _g.iroot(12,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: n must be > 0
>>> _g.iroot(12,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: n must be > 0
>>> _g.iroot(-12,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: iroot() of negative number
>>> _g.digits(3.14)
('31400000000000001', 1, 53)
>>> _g.digits(3,'peep')
Traceback (most recent call last):
  ...
TypeError: digits() requires 'int' argument for base
>>> _g.digits(3.14,'peep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.digits(3,'peep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: digits() requires 'int' argument for base
>>> _g.mpz(3).digits('bu')
Traceback (most recent call last):
  ...
TypeError: digits() requires 'int' argument for base
>>> _g.mpfr(3).digits('bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.qdiv(3,_g.mpq(1))
mpz(3)
>>> _g.qdiv(3,_g.mpz(1))
mpz(3)
>>> _g.qdiv(3,_g.mpfr(1))
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
TypeError: second argument cannot be converted to 'mpq'
>>> _g.qdiv(3,4,5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes at most 2 arguments (3 given)
>>> _g.qdiv('bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: first argument cannot be converted to 'mpq'
>>> _g.qdiv(1.0,1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: first argument cannot be converted to 'mpq'
>>> _g.qdiv(1,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: division or modulo by zero in qdiv
>>> _g.f2q(-1.0)
mpz(-1)
>>> _g.mpz(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes at most 2 arguments (3 given)
>>> _g.mpz('bi','bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.mpz('bi',99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for mpz() must be 0 or in the interval 2 ... 62
>>> _g.mpz(1,2)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpz() with non-string argument needs exactly 1 argument
>>> _g.mpz(None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpz() requires numeric or string argument
>>> _g.mpq(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpq() requires 0, 1 or 2 arguments
>>> _g.mpq('bi','bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.mpq('bi',99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for mpq() must be 0 or in the interval 2 ... 62
>>> _g.mpq(None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpq() requires numeric or string argument
>>> _g.mpq(1,None)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpq() requires numeric or string argument
>>> _g.mpq(1,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: zero denominator in 'mpq'
>>> _g.mpfr()
mpfr('0.0')
>>> _g.mpfr(1,'bo')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.mpfr(1,-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: precision for mpfr() must be >= 0
>>> _g.mpfr('ba',0,'bu')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: an integer is required
>>> _g.mpfr('ba',0,99)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: base for mpfr() must be 0 or in the interval 2 ... 62
>>> _g.mpfr(1,2,3)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: function takes at most 2 arguments (3 given)
>>> +_g.mpz(1)
mpz(1)
>>> +_g.mpfr(1)
mpfr('1.0')
>>> +_g.mpq(1)
mpq(1,1)
>>> _g.mpz(2)**-2
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: pow() exponent cannot be negative
>>> _g.mpz(2)**_g.mpz(1000000000000000L*1000000000000000L)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: pow() outrageous exponent
>>> pow(_g.mpz(2),3,0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: pow() 3rd argument cannot be 0
>>> pow(_g.mpz(2),3,-5)
mpz(-2)
>>> pow(_g.mpq(2),3,-5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: mpq.pow() no modulo allowed
>>> a=10000000000L**2
>>> _g.mpq(2)**a
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: mpq.pow() outrageous exponent
>>> _g.mpq(2)**_g.mpq(1,a)
mpfr('1.0')
>>> _g.mpq(2)**0
mpq(1,1)
>>> _g.mpq(2)**-1
mpq(1,2)
>>> _g.mpq(2)**_g.mpq(1,2)
mpfr('1.4142135623730951')
>>> _g.mpq(-2)**_g.mpq(1,2)
mpfr('nan')
>>> _g.mpq(0)**_g.mpq(1,2)
mpfr('0.0')
>>> _g.mpq(0)**-1
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: mpq.pow() 0 base to negative exponent
>>> _g.mpq(-1)**-1
mpq(-1,1)
>>> _g.mpfr(9,100)**2
mpfr('81.0')
>>> _g.mpfr(9,100)**0.5
mpfr('3.0')
>>> _g.mpfr(9,100)**_g.mpfr(0.5)
mpfr('3.0')
>>> _g.mpfr(0)**2
mpfr('0.0')
>>> pow(_g.mpfr(2),3,-5)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: pow() 3rd argument not allowed unless all arguments are integers
>>> _g.mpz(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpz' and 'str'
>>> _g.mpz(1)+_g.mpfr(1)
mpfr('2.0')
>>> _g.mpz(1)+_g.mpq(1)
mpq(2,1)
>>> _g.mpq(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpq' and 'str'
>>> _g.mpfr(1)+'bu'
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: unsupported operand type(s) for +: 'mpfr' and 'str'
>>> _g.mpfr(1)+_g.mpq(2)
mpfr('3.0')
>>> divmod(_g.mpz(3),0)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ZeroDivisionError: division or modulo by zero
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
ValueError: fac() of negative number
>>> _g.fib(-1)
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
ValueError: Fibonacci of negative number
>>> _g.comb(7,-3)
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
>>> _g.jacobi(23, -34)
Traceback (most recent call last):
  ...
ValueError: y must be odd and >0
>>> _g.legendre(23, -34)
Traceback (most recent call last):
  ...
ValueError: y must be odd and >0
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
        print "Unit tests for gmpy2 (extra cover)"
        print "    on Python %s" % sys.version
        print "Testing gmpy2 {0}".format(_g.version())
        print "  Mutliple-precision library:   {0}".format(_g.mp_version())
        print "  Floating-point library:       {0}".format(_g.mpfr_version())
        print "  Complex library:              {0}".format(_g.mpc_version())
        print "  Caching Values: (Number)      {0}".format(_g.get_cache()[0])
        print "  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1])

    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat:
        print
        print "Overall results for cvr:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)
