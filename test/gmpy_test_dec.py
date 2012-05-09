# partial unit test for gmpy/decimal interoperability
# note: broken in Python 2.4.0 due to a 2.4.0 bug, please update to 2.4.1
#       or better to allow decimal/most-anything-else interoperability!-)
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> filter(lambda x: not x.startswith('__'), dir(f))
['_copy', 'binary', 'ceil', 'digits', 'f2q', 'floor', 'getprec', 'getrprec', 'qdiv', 'reldiff', 'round', 'setprec', 'sign', 'sqrt', 'trunc']
>>>
'''
try: import decimal as _d
except ImportError: _d = None

import gmpy as _g, doctest, sys
__test__={}
f=_g.mpf('123.456')
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
>>> print _g.mpf(d)
12.34
>>> print f+d
135.796
>>> print d+f
135.796
>>> print q+d
801.463
>>> print d+q
801.463
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
>>> _g.getrprec(d)
53
>>> _g.fsqrt(d)==_g.mpf(d).sqrt()
1
>>>
'''

def _test(chat=None):
    python_version = sys.version_info[:3]
    if python_version == (2, 4, 0):
        print "You're using Python 2.4.0, which does not allow interoperability"
        print "  between decimal and other types (due to a bug fixed in 2.4.1)"
        print "  No point in testing, therefore -- please upgrade your Python!"
        return 0, 0
    if chat:
        print "Unit tests for gmpy 1.15 (decimal interoperation)"
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

    if not _d:
        if chat:
            print "Can't test, since can't import decimal"
        return 0, 0
    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat:
        print
        print "Overall results for dec:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)
