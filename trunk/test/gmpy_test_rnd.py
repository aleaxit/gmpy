# partial unit test for gmpy2 rand functionality
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> r
<built-in function rand>
>>>
'''

import gmpy2 as _g, doctest,sys
__test__={}
r = _g.rand

__test__['rand']=\
r'''
>>> r('error',1,2,3)
Traceback (most recent call last):
   ...
TypeError: function takes exactly 2 arguments (4 given)
>>> r('save')
Traceback (most recent call last):
   ...
RuntimeError: can't save before init
>>> r('init',20)
>>> r('unkn',99)
Traceback (most recent call last):
   ...
ValueError: unknown option 'unkn'
>>> r('seed',1234)
>>> for i in range(5):
...     print r('next',100),
...     if i==4: print
...
21 75 63 28 27
>>> alis=list("proktelnu")
>>> for i in range(10):
...     r('shuf',alis)
...     print ''.join(alis)
...
rtoulpnke
eoturlknp
plnuetokr
ekoprulnt
kpoutnrel
rutoneklp
ukeptnorl
onkrlpteu
lknteropu
enrkutlpo
>>> sav=r('save')
>>> print sav
774447212137
>>> for i in range(5):
...     r('shuf',alis)
...     print ''.join(alis)
...
elnuortpk
enutolpkr
eropulntk
plroutenk
ekonrtplu
>>> r('seed',sav)
>>> for i in range(5):
...     r('shuf',alis)
...     print ''.join(alis)
...
epkruotln
ekrtuplno
eoulrpktn
lpourtekn
enukotlpr
>>> r('seed',sav)
>>> for i in range(3):
...     print float(r('floa'))
...
0.44833278656
0.547296524048
0.895370483398
>>> r('seed',sav)
>>> for i in range(3):
...     print float(r('floa',6))
...
0.484375
0.90625
0.75
>>> r('seed',sav)
>>> for i in range(3):
...     print _g.f2q(r('floa',6),-6)
...
15/31
9/10
3/4
>>> r('seed',sav)
>>> for i in range(3):
...     print _g.f2q(r('floa',6))
...
31/64
29/32
3/4
>>> r('seed',sav)
>>> for i in range(5):
...     r('shuf',alis)
...     print ''.join(alis)
...
elnorutpk
enotrlpku
eurpolntk
plurotenk
ekrnutplo
>>> try: r('shuf','astring')
... except TypeError, e: print int("does not support item assignment" in str(e))
1
>>> r('shuf',23)
Traceback (most recent call last):
   ...
TypeError: 'shuf' needs mutable sequence
>>>
'''

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy2 (rnd functionality)"
['ModeMPFR', 'ModePython', 'RoundAwayZero', 'RoundDown', 'RoundToNearest', 'RoundToZero', 'RoundUp', 'acos', 'acosh', 'add', 'agm', 'ai', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 'binary', 'bincoef', 'bit_clear', 'bit_flip', 'bit_length', 'bit_mask', 'bit_scan0', 'bit_scan1', 'bit_set', 'bit_test', 'cbrt', 'cdiv', 'cdiv2exp', 'cdivmod', 'cdivmod2exp', 'ceil', 'clear_erangeflag', 'clear_flags', 'clear_inexactflag', 'clear_nanflag', 'clear_overflow', 'clear_underflow', 'cmod', 'cmod2exp', 'comb', 'const_catalan', 'const_euler', 'const_log2', 'const_pi', 'cos', 'cosh', 'cot', 'coth', 'csc', 'csch', 'denom', 'digamma', 'digits', 'div', 'divexact', 'divm', 'eint', 'erf', 'erfc', 'exp', 'exp10', 'exp2', 'expm1', 'f2q', 'fac', 'factorial', 'fdiv', 'fdiv2exp', 'fdivmod', 'fdivmod2exp', 'fib', 'fib2', 'floor', 'fma', 'fmod', 'fmod2exp', 'fms', 'gamma', 'gcd', 'gcdext', 'get_cache', 'get_emax', 'get_emax_max', 'get_emin', 'get_emin_min', 'get_max_precision', 'get_mode', 'get_precision', 'get_rounding', 'get_ternary', 'hamdist', 'hypot', 'inf', 'invert', 'is_erangeflag', 'is_even', 'is_inexactflag', 'is_inf', 'is_nan', 'is_nanflag', 'is_number', 'is_odd', 'is_overflow', 'is_power', 'is_prime', 'is_regular', 'is_square', 'is_underflow', 'is_zero', 'j0', 'j1', 'jacobi', 'jn', 'kronecker', 'lcm', 'legendre', 'lgamma', 'li2', 'license', 'lngamma', 'log', 'log10', 'log1p', 'log2', 'lucas', 'lucas2', 'max', 'min', 'mp_limbsize', 'mp_version', 'mpf', 'mpfr_version', 'mpq', 'mpz', 'mul', 'nan', 'next_above', 'next_below', 'next_prime', 'next_toward', 'numdigits', 'numer', 'pack', 'pi', 'popcount', 'pow', 'qdiv', 'rec_sqrt', 'reldiff', 'remove', 'root', 'rootrem', 'round', 'sec', 'sech', 'set_cache', 'set_debug', 'set_emax', 'set_emin', 'set_erangeflag', 'set_fcoform', 'set_inexactflag', 'set_mode', 'set_nanflag', 'set_overflow', 'set_precision', 'set_rounding', 'set_underflow', 'sign', 'sin', 'sin_cos', 'sinh', 'sinh_cosh', 'sqrt', 'sqrtrem', 'square', 'sub', 'tan', 'tanh', 'tdiv', 'tdiv2exp', 'tdivmod', 'tdivmod2exp', 'tmod', 'tmod2exp', 'trunc', 'unpack', 'version', 'xbit_mask', 'xmpz', 'y0', 'y1', 'yn', 'zero', 'zeta']

    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat:
        print
        print "Overall results for rnd:"
    return doctest.master.summarize(chat)


if __name__=='__main__':
    _test(1)

