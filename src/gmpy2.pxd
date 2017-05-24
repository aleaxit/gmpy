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

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPZ_New() declared below.
    cdef (MPZ_Object *) GMPy_MPZ_New(void *)

    # WARNING: this function might cause memory leak if not used
    # appropriately. Prefer MPQ_New() declared below.
    cdef (MPQ_Object *) GMPy_MPQ_New(void *)

    # TODO
    # cdef (MPFR_Object *) GMPy_MPFR_New(void *)

    # TODO
    # cdef (MPC_Object *) GMPy_MPC_New(void *)

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

# TODO
# Build a gmpy2 mpfr from a mpfr
# cdef inline GMPy_MPQ_From_mpfr(mpfr_srcptr x)

# TODO
# Build a gmpy2 mpc from a mpc
# cdef inline GMPy_MPC_From_mpc(mpc_srcptr z)

# TODO
# Build a gmpy2 mpc from a real part mpfr and an imaginary part mpfr
# cdef inline GMPy_MPC_From_mpfr(mpfr_srcptr re, mpfr_scrptr im)
