from cpython.object cimport PyObject

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


cdef extern from "mpc.h":
    # mpc complexes
    ctypedef struct __mpc_struct:
        pass
    ctypedef __mpc_struct mpc_t[1];
    ctypedef __mpc_struct *mpc_ptr;
    ctypedef const __mpc_struct *mpc_srcptr;


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

    # construct a new gmpy2 mpz
    # (equivalent to mpz.__new__(mpz)
    cdef (MPZ_Object *) GMPy_MPZ_New(void *)

    # construct a new gmpy2 mpq
    # (equivalent of mpq.__new__(mpq)
    cdef (MPQ_Object *) GMPy_MPQ_New(void *)

    # TODO
    # construct a new gmpy2 mpfr
    # (equivalent of mpfr.__new__(mpfr)
    # cdef (MPFR_Object *) GMPy_MPFR_New(void *)

    # TODO
    # construct a new gmpy mpc
    # (equivalent of mpc.__new__(mpc)
    # cdef (MPC_Object *) GMPy_MPC_New(void *)

    # access to the mpz_t field of a gmpy2 mpz
    cdef mpz_t MPZ(MPZ_Object *)

    # access to the mpq_t field of a gmpy2 mpq
    cdef mpq_t MPQ(MPQ_Object *)

    # access to the mpfr_t field of a gmpy2 mpfr
    cdef mpfr_t MPFR(MPFR_Object *)

    # access to the mpc_t field of a gmpy2 mpc
    cdef mpc_t MPC(MPC_Object *)

# Build a gmpy2 mpz from a gmp mpz
cdef inline GMPy_MPZ_From_mpz(mpz_srcptr z):
    cdef MPZ_Object * res = GMPy_MPZ_New(NULL)
    mpz_set(MPZ(res), z)
    return <object> res

# Build a gmpy2 mpq from a gmp mpq
cdef inline GMPy_MPQ_From_mpq(mpq_srcptr q):
    cdef MPQ_Object * res = GMPy_MPQ_New(NULL)
    mpq_set(MPQ(res), q)
    return <object> res

# Build a gmpy2 mpq from gmp mpz numerator and denominator
cdef inline GMPy_MPQ_From_mpz(mpz_srcptr num, mpz_srcptr den):
    cdef MPQ_Object * res = GMPy_MPQ_New(NULL)
    mpq_set_num(MPQ(res), num)
    mpq_set_den(MPQ(res), den)
    return <object> res

# TODO
# Build a gmpy2 mpfr from a mpfr
# cdef inline GMPy_MPQ_From_mpfr(mpfr_srcptr x)

# TODO
# Build a gmpy2 mpc from a mpc
# cdef inline GMPy_MPC_From_mpc(mpc_srcptr z)

# TODO
# Build a gmpy2 mpc from a real part mpfr and an imaginary part mpfr
# cdef inline GMPy_MPC_From_mpfr(mpfr_srcptr re, mpfr_scrptr im)
