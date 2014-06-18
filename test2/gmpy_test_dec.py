# partial unit test for gmpy2/decimal interoperability
# relies on Tim Peters' "doctest.py" test-driver
r'''
>>> filter(lambda x: not x.startswith('__'), dir(f))
['as_integer_ratio', 'as_mantissa_exp', 'as_simple_fraction', 'conjugate', 'digits', 'imag', 'is_integer', 'precision', 'rc', 'real']
>>>
'''
try: import decimal as _d
except ImportError: _d = None

import gmpy2 as _g, doctest, sys
__test__={}
f=_g.mpfr('123.456')
q=_g.mpq('789123/1000')
z=_g.mpz('234')
if _d: d=_d.Decimal('12.34')

__test__['elemop']=\
r'''
>>> print _g.mpz(23) == _d.Decimal(23)
True
>>> print _g.mpz(d)
12
>>> print _g.mpq(d)
617/50
>>> print _g.mpfr(d)
12.34
>>> print f+d
135.79599999999999
>>> print d+f
135.79599999999999
>>> print q+d
801.46300000000008
>>> print d+q
801.46300000000008
>>> print z+d
246.34
>>> print d+z
246.34
>>> print _g.ceil(d)
13.0
>>> print _g.floor(d)
12.0
>>> print _g.trunc(d)
12.0
>>> _g.mpfr(d).precision
53
>>>
'''

def _test(chat=None):
    python_version = sys.version_info[:3]
    if chat:
        print "Unit tests for gmpy2 (decimal interoperation)"
        print "    on Python %s" % sys.version
        print "Testing gmpy2 {0}".format(_g.version())
        print "  Mutliple-precision library:   {0}".format(_g.mp_version())
        print "  Floating-point library:       {0}".format(_g.mpfr_version())
        print "  Complex library:              {0}".format(_g.mpc_version())
        print "  Caching Values: (Number)      {0}".format(_g.get_cache()[0])
        print "  Caching Values: (Size, limbs) {0}".format(_g.get_cache()[1])

    if not _d:
        if chat:
            print "Can't test, since can't import decimal"
        return 0, 0
    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0, optionflags=doctest.IGNORE_EXCEPTION_DETAIL)

    if chat:
        print
        print "Overall results for dec:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)
