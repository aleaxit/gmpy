from cpython.object cimport PyObject

cdef extern from "Python.h":
    void Py_DECREF(PyObject *)

cdef extern from "gmp.h":
    # gmp integers
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

    # gmp rationnals
    ctypedef struct __mpq_struct:
        pass
    ctypedef __mpq_struct mpq_t[1]
    ctypedef __mpq_struct *mpq_ptr
    ctypedef const __mpq_struct *mpq_srcptr

    void mpz_set(mpz_t rop, mpz_t op)
    void mpq_set(mpq_ptr rop, mpq_srcptr op)
    void mpq_set_num (mpq_t rational, mpz_t numerator)
    void mpq_set_den (mpq_t rational, mpz_t denominator)


cdef extern from "mpfr.h":
    # mpfr reals
    ctypedef struct __mpfr_struct:
        pass
    ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct *mpfr_ptr
    ctypedef const __mpfr_struct *mpfr_srcptr

    ctypedef long mpfr_prec_t

    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN
        MPFR_RNDZ
        MPFR_RNDU
        MPFR_RNDD
        MPFR_RNDA
        MPFR_RNDF
        MPFR_RNDNA
        GMP_RNDN
        GMP_RNDZ
        GMP_RNDU
        GMP_RNDD

    mpfr_prec_t mpfr_get_prec (mpfr_t x)
    int mpfr_set (mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd)

cdef extern from "mpc.h":
    # mpc complexes
    ctypedef struct __mpc_struct:
        pass
    ctypedef __mpc_struct mpc_t[1];
    ctypedef __mpc_struct *mpc_ptr;
    ctypedef const __mpc_struct *mpc_srcptr;
    ctypedef enum mpc_rnd_t:
        MPC_RNDNN
        MPC_RNDNZ
        MPC_RNDNU
        MPC_RNDND
        MPC_RNDZN
        MPC_RNDZZ
        MPC_RNDZU
        MPC_RNDZD
        MPC_RNDUN
        MPC_RNDUZ
        MPC_RNDUU
        MPC_RNDUD
        MPC_RNDDN
        MPC_RNDDZ
        MPC_RNDDU
        MPC_RNDDD

    mpfr_prec_t mpc_get_prec (mpc_srcptr x)
    void mpc_get_prec2 (mpfr_prec_t *pr, mpfr_prec_t *pi, mpc_srcptr x)
    int  mpc_set (mpc_ptr rop, mpc_srcptr op, mpc_rnd_t rnd)
    int  mpc_set_fr_fr (mpc_ptr rop, mpfr_srcptr rp, mpfr_srcptr ip, mpc_rnd_t rnd)

cdef extern from "gmpy2/gmpy2.h":
    # initialize the C-API
    # need to be called before any other functions
    cdef int import_gmpy2()

    # object types
    ctypedef struct MPZ_Object:
        pass
    ctypedef struct MPQ_Object:
        pass
    ctypedef struct MPFR_Object:
        pass
    ctypedef struct MPC_Object:
        pass

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPZ_New() declared below.
    cdef (MPZ_Object *) GMPy_MPZ_New(void *)

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPQ_New() declared below.
    cdef (MPQ_Object *) GMPy_MPQ_New(void *)

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPFR_New() declared below.
    cdef (MPFR_Object *) GMPy_MPFR_New(mpfr_prec_t prec, void *)

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPC_New() declared below.
    cdef (MPC_Object *) GMPy_MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec, void *)

    # access to the mpz_t field of a gmpy2 mpz
    cdef mpz_t MPZ(MPZ_Object *)

    # access to the mpq_t field of a gmpy2 mpq
    cdef mpq_t MPQ(MPQ_Object *)

    # access to the mpfr_t field of a gmpy2 mpfr
    cdef mpfr_t MPFR(MPFR_Object *)

    # access to the mpc_t field of a gmpy2 mpc
    cdef mpc_t MPC(MPC_Object *)

    # check if "param" is a MPZ_Object
    cdef bint MPZ_Check(PyObject *)

    # check if "param" is a MPQ_Object
    cdef bint MPQ_Check(PyObject *)

    # check if "param" is a MPFR_Object
    cdef bint MPFR_Check(PyObject *)

    # check if "param" is a MPC_Object
    cdef bint MPC_Check(PyObject *)

# Return a new gmpy2 mpz object
# equivalent to mpz.__new__(mpz)
cdef inline MPZ_New():
    res = <object> GMPy_MPZ_New(NULL)  # Cython increases refcount!
    Py_DECREF(<PyObject *> res)
    return res

# Return a new gmpy2 mpq object
# equivalent to mpq.__new__(mpq)
cdef inline MPQ_New():
    res = <object> GMPy_MPQ_New(NULL)  # Cython increases refcount!
    Py_DECREF(<PyObject *> res)
    return res

# Return a new gmpy2 mpfr object
# equivalent to mpfr.__new__(mpfr)
cdef inline MPFR_New(mpfr_prec_t prec):
    res = <object> GMPy_MPFR_New(prec, NULL)  # Cython increases refcount!
    Py_DECREF(<PyObject *> res)
    return res

# Return a new gmpy2 mpc object
# equivalent to mpc.__new__(mpc)
cdef inline MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec):
    res = <object> GMPy_MPC_New(rprec, iprec, NULL)  # Cython increases refcount!
    Py_DECREF(<PyObject *> res)
    return res

# Build a gmpy2 mpz from a gmp mpz
cdef inline GMPy_MPZ_From_mpz(mpz_srcptr z):
    res = MPZ_New()
    mpz_set(MPZ(<MPZ_Object *> res), z)
    return res

# Build a gmpy2 mpq from a gmp mpq
cdef inline GMPy_MPQ_From_mpq(mpq_srcptr q):
    res = MPQ_New()
    mpq_set(MPQ(<MPQ_Object *> res), q)
    return res

# Build a gmpy2 mpq from gmp mpz numerator and denominator
cdef inline GMPy_MPQ_From_mpz(mpz_srcptr num, mpz_srcptr den):
    res = MPQ_New()
    mpq_set_num(MPQ(<MPQ_Object *> res), num)
    mpq_set_den(MPQ(<MPQ_Object *> res), den)
    return res

# Build a gmpy2 mpfr from a mpfr
cdef inline GMPy_MPFR_From_mpfr(mpfr_srcptr x):
    res = MPFR_New(mpfr_get_prec(x))
    mpfr_set(MPFR(<MPFR_Object *> res), x, MPFR_RNDN)
    return res

# Build a gmpy2 mpc from a mpc
cdef inline GMPy_MPC_From_mpc(mpc_srcptr c):
    cdef mpfr_prec_t pr
    cdef mpfr_prec_t pi
    mpc_get_prec2(&pr, &pi, c)
    res = MPC_New(pr, pi)
    mpc_set(MPC(<MPC_Object *> res), c, MPC_RNDNN)
    return res

# Build a gmpy2 mpc from a real part mpfr and an imaginary part mpfr
cdef inline GMPy_MPC_From_mpfr(mpfr_srcptr re, mpfr_srcptr im):
    res = MPC_New(mpfr_get_prec(re), mpfr_get_prec(im))
    mpc_set_fr_fr (MPC(<MPC_Object *> res), re, im, MPC_RNDNN)
    return res
