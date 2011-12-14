# partial unit test for gmpy mpq functionality
# relies on Tim Peters' "doctest.py" test-driver

r'''
>>> dir(a)
['__abs__', '__add__', '__bool__', '__class__', '__delattr__', '__divmod__', '__doc__', '__eq__', '__float__', '__floordiv__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__int__', '__le__', '__lt__', '__mod__', '__mul__', '__ne__', '__neg__', '__new__', '__pos__', '__pow__', '__radd__', '__rdivmod__', '__reduce__', '__reduce_ex__', '__repr__', '__rfloordiv__', '__rmod__', '__rmul__', '__rpow__', '__rsub__', '__rtruediv__', '__setattr__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '_copy', 'binary', 'denom', 'denominator', 'digits', 'numer', 'numerator', 'qdiv', 'sign']
>>>
'''

import gmpy as _g, doctest,sys
import fractions
F=fractions.Fraction

__test__={}
a=_g.mpq('123/456')
b=_g.mpq('789/123')
af=F(123,456)
bf=F(789,123)

__test__['compat']=\
r'''
>>> a==af
True
>>> af==a
True
>>> a < af
False
>>> a <= af
True
>>> a > af
False
>>> a >= af
True
>>> af < a
False
>>> af <= a
True
>>> af > a
False
>>> af >= a
True
>>> a+bf
mpq(41657,6232)
>>> divmod(123*a, b) == divmod(123*af, bf)
True
>>> divmod(-23*a, b) == divmod(-23*af, bf)
True
>>> divmod(a+17, b-23) == divmod(af+17, bf-23)
True
>>> divmod(-a, -b) == divmod(-af, -bf)
True
'''

__test__['elemop']=\
r'''
>>> a+b
mpq(41657,6232)
>>> a-b
mpq(-38295,6232)
>>> a*b
mpq(263,152)
>>> a/b
mpq(1681,39976)
>>> a//b
mpz(0)
>>> a//-b
mpz(-1)
>>> -a//b
mpz(-1)
>>> -a//-b
mpz(0)
>>> b+a
mpq(41657,6232)
>>> b-a
mpq(38295,6232)
>>> b*a
mpq(263,152)
>>> b/a
mpq(39976,1681)
>>> b//a
mpz(23)
>>> a+1
mpq(193,152)
>>> 1+a
mpq(193,152)
>>> a-1
mpq(-111,152)
>>> 1-a
mpq(111,152)
>>> a*1
mpq(41,152)
>>> 1*a
mpq(41,152)
>>> a/1
mpq(41,152)
>>> 1/a
mpq(152,41)
>>> a % b
mpq(41,152)
>>> a % -b
mpq(-38295,6232)
>>> 2*a % 7*b
mpq(263,76)
>>> -a
mpq(-41,152)
>>> abs(-a)
mpq(41,152)
>>> _g.qsign(b-a)
1
>>> _g.qsign(b-b)
0
>>> _g.qsign(a-b)
-1
>>> a.sign()
1
>>> (-a).sign()
-1
>>> z=b-b; z.sign()
0
>>> a.numer() == a.numerator
True
>>> a.denom() == a.denominator
True
>>> an=a.numer(); ad=a.denom();
>>> an==0 or 1==a*_g.mpq(ad,an)
1
>>> bn=b.numer(); bd=b.denom();
>>> bn==0 or 1==b*_g.mpq(bd,bn)
1
>>> zn=z.numer(); zd=z.denom();
>>> zn==0 or 1==z*_g.mpq(zd,zn)
1
>>> (a+b) == _g.mpq(an*bd+ad*bn,ad*bd)
1
>>> (a+z) == _g.mpq(an*zd+ad*zn,ad*zd)
1
>>> (a+a) == _g.mpq(an*ad+ad*an,ad*ad)
1
>>> import pickle
>>> pickle.loads(pickle.dumps(_g.mpq(1234,6789)))
mpq(1234,6789)
>>>
'''

__test__['cmpr']=\
r'''
>>> c=_g.mpq(a)
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
>>> a>b
0
>>> a<b
1
>>> not _g.mpq(0)
1
>>> not a
0
>>> a>1
0
>>> a>1.0
0
>>> a<1
1
>>> a<1.0
1
>>> a==1
0
>>> a==1.0
0
>>> int(1/a)
3
>>> int(-1/a)
-3
>>>
'''

__test__['format']=\
r'''
>>> str(a)
'41/152'
>>> repr(a)
'mpq(41,152)'
>>> a==eval(repr(a),_g.__dict__)
1
>>> str(-a)
'-41/152'
>>> repr(-a)
'mpq(-41,152)'
>>> (-a)==eval(repr(-a),_g.__dict__)
1
>>> _g.set_tagoff(0)
1
>>> a
gmpy.mpq(41,152)
>>> _g.mpq('12.34')
gmpy.mpq(617,50)
>>> _g.set_tagoff(1)
0
>>> for i in range(1,7):
...    for j in range(3,10):
...       if _g.mpq(i,j) != _g.mpq("%d/%d"%(i,j)):
...          print('er1:',i,j); break
...       aa=_g.mpq(i,j); ai=aa.numer(); aj=aa.denom()
...       if aj!=1 and str(aa) != ("%d/%d"%(ai,aj)):
...          print('er2:',i,j,str(aa),("%d/%d"%(ai,aj))); break
...       if aj==1 and str(aa) != ("%d"%ai):
...          print('er3:',i,j,str(aa),"%d"%ai); break
...       if aj!=1 and repr(aa) != ("mpq(%d,%d)"%(ai,aj)):
...          print('er4:',i,j,repr(aa),("mpq(%d,%d)"%(ai,aj))); break
...       if aj==1 and repr(aa) != ("mpq(%d,%d)"%(ai,aj)):
...          print('er5:',i,j,repr(aa),"mpq(%d,%d)"%(ai,aj)); break
>>> fmo='_g.mpq('+hex(a.numer())+','+hex(a.denom())+')'
>>> fmo
'_g.mpq(0x29,0x98)'
>>> eval(fmo)==a
1
>>> fmo='_g.mpq("'+a.numer().digits(30)+'/'+a.denom().digits(30)+'",30)'
>>> fmo
'_g.mpq("1b/52",30)'
>>> eval(fmo)==a
1
>>> _g.qdigits(a,30)
'1b/52'
>>> a.digits(30)
'1b/52'
>>> _g.mpq(1000*1000*1000*1000*1000*1000*1000,23)
mpq(1000000000000000000000,23)
>>> _g.mpq(23,1000*1000*1000*1000*1000*1000*1000)
mpq(23,1000000000000000000000)
>>> _g.mpq(23**15,1000**7)
mpq(266635235464391245607,1000000000000000000000)
>>> _g.qbinary('pep')
Traceback (most recent call last):
  File "<stdin>", line 1, in ?
TypeError: argument can not be converted to mpq
>>> x=_g.mpq('234/567')
>>> del x
>>> _g.mpq('7788')
mpq(7788,1)
>>> _g.mpq('12.34')
mpq(617,50)
'''

__test__['binio']=\
r'''
>>> ba=a.binary()
>>> len(ba)
6
>>> for i in range(len(ba)):
...     print(ba[i])
...
1
0
0
0
41
152
>>> _g.mpq(ba,256)==a
1
>>> ba == _g.qbinary(a)
1
>>> ba=(-a).binary()
>>> len(ba)
6
>>> for i in range(len(ba)):
...     print(ba[i])
...
1
0
0
128
41
152
>>> _g.mpq(ba,256)==-a
1
>>>
'''

__test__['power']=\
r'''
>>> _g.mpq(2,3)**3
mpq(8,27)
>>> _g.mpq(8,27)**_g.mpq('2/3')
mpq(4,9)
>>> _g.mpq(2,3)**-3
mpq(27,8)
>>> _g.mpq(8,27)**_g.mpq('-2/3')
mpq(9,4)
>>> print(float("%.14f" % _g.mpf('0.2')**2))
0.04
>>> print(float(_g.mpf('0.2')**-2))
25.0
>>> _g.mpq(3)**3 == _g.mpz(3)**3
True
>>> (a**-7) == 1/(a**7)
True
>>> (b**5) == 1/(b**-5)
True
>>>
'''

__test__['qdiv']=\
r'''
>>> _g.qdiv(12,2)
mpz(6)
>>> _g.qdiv(12,5)
mpq(12,5)
>>> a is a.qdiv()
1
>>> a is a.qdiv(1)
1
>>> a is a.qdiv(2)
0
>>> x=a.numer()
>>> x is _g.qdiv(x)
1
>>> x is _g.qdiv(x,1)
1
>>> x is _g.qdiv(x,2)
0
>>> y=_g.mpq(4,1)
>>> y is y.qdiv()
0
>>> y == y.qdiv()
1
>>>
'''

def _test(chat=None):
    if chat:
        print("Unit tests for gmpy 1.15 (mpq functionality)")
        print("    running on Python",sys.version)
        print()
        if _g.gmp_version():
            print("Testing gmpy %s (GMP %s), default caching (%s, %s)" % (
                (_g.version(), _g.gmp_version(), _g.get_cache()[0],
                _g.get_cache()[1])))
        else:
            print("Testing gmpy %s (MPIR %s), default caching (%s, %s)" % (
                (_g.version(), _g.mpir_version(), _g.get_cache()[0],
                _g.get_cache()[1])))

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
        doctest.testmod(thismod, report=0)
    finally:
        sys.stdout = sav

    if chat:
        print()
        print("Overall results for mpq:")
    return doctest.master.summarize(chat)


if __name__=='__main__':
    _test(1)

