# partial unit test for gmpy's "set" functionality
# relies on Tim Peters' "doctest.py" test-driver
# test-version 1.02
r'''
>>> dir(s)
['add', 'getimax']
>>>
'''
import gmpy as _g, doctest, sys
__test__={}
s=_g.intset(52)

__test__['elemop']=\
r'''
>>> type(s)
<type 'intset'>
>>> s.getimax()
52
>>> len(s)
0
>>> s.add(23)
>>> s.add('x')
Traceback (most recent call last):
  ...
TypeError: an integer is required
>>> s.add(-3)
Traceback (most recent call last):
  ...
ValueError: n must be >= 0
>>> s.add(52)
Traceback (most recent call last):
  ...
ValueError: n must be < the intset's getimax()
>>> s.add(52431)
Traceback (most recent call last):
  ...
ValueError: n must be < the intset's getimax()
>>> len(s)
1
>>> s.add(45)
>>> len(s)
2
>>> s.add(23)
>>> len(s)
2
>>> s.add(45)
>>> len(s)
2
>>> 23 in s
True
>>> 32 in s
False
>>> 45 in s
True
>>> 54 in s
False
>>> -33 in s
False
>>> 'x' in s
Traceback (most recent call last):
  ...
TypeError: an integer is required
>>> (sys.maxint+333) in s
False
>>> (-sys.maxint-333) in s
False
'''

def _test(chat=None):
    if chat:
        print "Unit tests for gmpy 1.03 release candidate (set functionality)"
        print "    running on Python", sys.version
        print
        print "Testing gmpy %s (GMP %s)" % (_g.version(), _g.gmp_version())
    thismod = sys.modules.get(__name__)
    doctest.testmod(thismod, report=0)

    if chat:
        print
        print "Overall results for set:"
    return doctest.master.summarize(chat)

if __name__=='__main__':
    _test(1)

