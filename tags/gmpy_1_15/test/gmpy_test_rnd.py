# partial unit test for gmpy rand functionality
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> r
<built-in function rand>
>>>
'''

import gmpy as _g, doctest,sys
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

# adapt to python 2.3's slightly different error message in an exception
import sys
if sys.version<'2.4':
    __test__['rand'] = __test__['rand'].replace("does not", "doesn't")

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy 1.15 (rnd functionality)"
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

    if chat:
        print
        print "Overall results for rnd:"
    return doctest.master.summarize(chat)


if __name__=='__main__':
    _test(1)

