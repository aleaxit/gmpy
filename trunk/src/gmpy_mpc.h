/*
  gmpy_mpc C API extension header file.

  Provide interface to the MPC (Multiple Precision Complex) library.

  Version 2.00, April 2011 (created) casevh

  This file is expected to be included from gmpy.h
 */

#if defined(MS_WIN32) && defined(_MSC_VER)
#  pragma comment(lib,"mpc.lib")
#endif

typedef struct {
    mpob ob;
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
    int round_mode;
} PympcObject;

static PyTypeObject Pympc_Type;
#define Pympc_AS_MPC(obj) (((PympcObject *)(obj))->c)
#define Pympc_Check(v) (((PyObject*)v)->ob_type == &Pympc_Type)
