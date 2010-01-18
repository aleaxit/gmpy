/*
  gmpy C API extension header file.
  Part of Python's gmpy module since version 0.4

  Created by Pearu Peterson <pearu@cens.ioc.ee>, November 2000.
  Edited by A. Martelli <aleaxit@yahoo.com>, December 2000.
  Edited by Case Van Horsen <casevh@gmail.com>, 2009, 2010.

  Version 1.02, February 2007.
  Version 1.03, June 2008
  Version 1.04, June 2008 (no changes)
  Version 1.05, February 2009 (support MPIR)
  Version 1.20, January 2010 (remove obsolete MS hacks, added
                PyInt_FromSize_t) casevh
 */

#ifndef Py_GMPYMODULE_H
#define Py_GMPYMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
   /* so one won't need to link explicitly to gmp.lib...: */
#  if defined(MPIR)
#    pragma comment(lib,"mpir.lib")
#  else
#    pragma comment(lib,"gmp.lib")
#  endif
#  define isnan _isnan
#  define isinf !_finite
#  define USE_ALLOCA 1
#  define inline __inline
#endif

#if defined MPIR
#  include "mpir.h"
#else
#  include "gmp.h"
#endif

#ifndef Py_TPFLAGS_HAVE_INDEX
#  define Py_TPFLAGS_HAVE_INDEX 0
#endif

#ifdef __GNUC__
#define USE_ALLOCA 1
#endif

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#   ifdef _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    if HAVE_ALLOCA_H
#     include <alloca.h>
#    else
       char *alloca ();
#    endif
#   endif
# endif
#endif

#define ALLOC_THRESHOLD 8192

#ifdef USE_ALLOCA
#define TEMP_ALLOC(B, S) \
    if(S < ALLOC_THRESHOLD) { \
        B = alloca(S); \
    } else { \
        if(!(B = PyMem_Malloc(S))) { \
            PyErr_NoMemory(); \
            return NULL; \
        } \
    }
#define TEMP_FREE(B, S) if(S >= ALLOC_THRESHOLD) PyMem_Free(B)
#else
#define TEMP_ALLOC(B, S) \
    if(!(B = PyMem_Malloc(S)))  { \
        PyErr_NoMemory(); \
        return NULL; \
    }
#define TEMP_FREE(B, S) PyMem_Free(B)
#endif

/* Various defs to mask differences between Python versions. */

#define Py_RETURN_NOTIMPLEMENTED\
    return Py_INCREF(Py_NotImplemented), Py_NotImplemented

#ifndef Py_SIZE
#define Py_SIZE(ob)     (((PyVarObject*)(ob))->ob_size)
#endif

#ifndef Py_TYPE
#define Py_TYPE(ob)     (((PyObject*)(ob))->ob_type)
#endif

#ifndef Py_REFCNT
#define Py_REFCNT(ob)   (((PyObject*)(ob))->ob_refcnt)
#endif

#if PY_VERSION_HEX < 0x02050000
#  error "GMPY 1.2 requires Python 2.5 or later."
#endif

/* Header file for gmpy */
typedef struct {
    PyObject_HEAD
    /* PyObject* callable; */
    /* long flags; */
} mpob;
typedef struct {
    mpob ob;
    mpz_t z;
} PympzObject;
typedef struct {
    mpob ob;
    mpq_t q;
} PympqObject;
typedef struct {
    mpob ob;
    mpf_t f;
    unsigned int rebits;
} PympfObject;

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)
#define Pympq_AS_MPQ(obj) (((PympqObject *)(obj))->q)
#define Pympf_AS_MPF(obj) (((PympfObject *)(obj))->f)

/* This section is used when compiling gmpy.c */
static PyTypeObject Pympz_Type;
#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)
static PyTypeObject Pympq_Type;
#define Pympq_Check(v) (((PyObject*)v)->ob_type == &Pympq_Type)
static PyTypeObject Pympf_Type;
#define Pympf_Check(v) (((PyObject*)v)->ob_type == &Pympf_Type)

//~ static PympzObject * Pympz_new(void);
//~ static void Pympz_dealloc(PympzObject *self);
//~ static int Pympz_convert_arg(PyObject *arg, PyObject **ptr);

//~ static PympqObject * Pympq_new(void);
//~ static void Pympq_dealloc(PympqObject *self);
//~ static int Pympq_convert_arg(PyObject *arg, PyObject **ptr);

//~ static PympfObject * Pympf_new(unsigned long bits);
//~ static void Pympf_dealloc(PympfObject *self);
//~ static int Pympf_convert_arg(PyObject *arg, PyObject **ptr);

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
