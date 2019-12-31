import gmpy2

class Z:
    def __mpz__(self): return gmpy2.mpz(2)
class Q:
    def __mpz__(self): return gmpy2.mpz(1)
    def __mpq__(self): return gmpy2.mpq(3,2)
class R:
    def __mpfr__(self): return gmpy2.mpfr(1.5)
class Cx:
    def __mpfr__(self): return gmpy2.mpfr(1.5)
    def __mpc__(self): return gmpy2.mpc(42,67)

z = Z()
q = Q()
r = R()
cx = Cx()

class A:
    def __mpz__(self): return gmpy2.mpz(1)
    def __mpq__(self): return gmpy2.mpq(3,2)
    def __mpfr__(self): return gmpy2.mpfr(1.5)
    def __mpc__(self): return gmpy2.mpc(42,67)
class B:
    def __mpz__(self): return 'hello'
    def __mpq__(self): return 'hello'
    def __mpfr__(self): return 'hello'
    def __mpc__(self): return 'hello'
class C:
    pass
class D:
    def __mpz__(self): raise TypeError
    def __mpq__(self): raise TypeError
    def __mpfr__(self): raise TypeError
    def __mpc__(self): raise TypeError

a = A()
b = B()
c = C()
d = D()
