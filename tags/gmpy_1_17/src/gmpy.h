/*
  gmpy C API extension header file.
  Part of Python's gmpy module since version 0.4

  Created by Pearu Peterson <pearu@cens.ioc.ee>, November 2000.
  Edited by A. Martelli <aleaxit@yahoo.com>, December 2000.
  Version 1.02, February 2007.
  Version 1.03, June 2008
  Version 1.04, June 2008 (no changes)
  Version 1.05, February 2009 (support MPIR)
 */

#ifndef Py_GMPYMODULE_H
#define Py_GMPYMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
/* the __MPN determination in stock gmp.h doesn't work, so...: */
#    define __MPN(x) __gmpn_##x
#    define _PROTO(x) x
#define inline __inline
#endif

#if defined MPIR
#include "mpir.h"
#else
#include "gmp.h"
#endif

/* ensure 2.5 compatibility */
#if PY_VERSION_HEX < 0x02050000
typedef int Py_ssize_t;
#define PY_FORMAT_SIZE_T ""
#endif
#ifndef Py_TPFLAGS_HAVE_INDEX
#define Py_TPFLAGS_HAVE_INDEX 0
#endif

#if PY_VERSION_HEX < 0x030200A3
typedef long Py_hash_t;
typedef unsigned long Py_uhash_t;
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
    size_t rebits;
} PympfObject;

/* #define MPOBCAL(obj) ((mpob*)obj)->callable */
/* #define MPOBFLA(obj) ((mpob*)obj)->flags */

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)
#define Pympq_AS_MPQ(obj) (((PympqObject *)(obj))->q)
#define Pympf_AS_MPF(obj) (((PympfObject *)(obj))->f)


#define Pympz_Type_NUM 0
#define Pympq_Type_NUM 1
#define Pympf_Type_NUM 2

  /* C API functions */

#define Pympz_new_NUM 3
#define Pympz_new_RETURN PympzObject *
#define Pympz_new_PROTO (void)

#define Pympq_new_NUM 4
#define Pympq_new_RETURN PympqObject *
#define Pympq_new_PROTO (void)

#define Pympf_new_NUM 5
#define Pympf_new_RETURN PympfObject *
#define Pympf_new_PROTO (size_t bits)

#define Pympz_dealloc_NUM 6
#define Pympz_dealloc_RETURN void
#define Pympz_dealloc_PROTO (PympzObject *self)

#define Pympq_dealloc_NUM 7
#define Pympq_dealloc_RETURN void
#define Pympq_dealloc_PROTO (PympqObject *self)

#define Pympf_dealloc_NUM 8
#define Pympf_dealloc_RETURN void
#define Pympf_dealloc_PROTO (PympfObject *self)

#define Pympz_convert_arg_NUM 9
#define Pympz_convert_arg_RETURN int
#define Pympz_convert_arg_PROTO (PyObject *arg, PyObject **ptr)

#define Pympq_convert_arg_NUM 10
#define Pympq_convert_arg_RETURN int
#define Pympq_convert_arg_PROTO (PyObject *arg, PyObject **ptr)

#define Pympf_convert_arg_NUM 11
#define Pympf_convert_arg_RETURN int
#define Pympf_convert_arg_PROTO (PyObject *arg, PyObject **ptr)

/* Total number of C API pointers */
#define Pygmpy_API_pointers 12

#ifdef GMPY_MODULE
/* This section is used when compiling gmpy.c */
static PyTypeObject Pympz_Type;
#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)
static PyTypeObject Pympq_Type;
#define Pympq_Check(v) (((PyObject*)v)->ob_type == &Pympq_Type)
static PyTypeObject Pympf_Type;
#define Pympf_Check(v) (((PyObject*)v)->ob_type == &Pympf_Type)

static Pympz_new_RETURN Pympz_new Pympz_new_PROTO;
static Pympz_dealloc_RETURN Pympz_dealloc Pympz_dealloc_PROTO;
static Pympz_convert_arg_RETURN Pympz_convert_arg Pympz_convert_arg_PROTO;
static Pympq_new_RETURN Pympq_new Pympq_new_PROTO;
static Pympq_dealloc_RETURN Pympq_dealloc Pympq_dealloc_PROTO;
static Pympq_convert_arg_RETURN Pympq_convert_arg Pympq_convert_arg_PROTO;
static Pympf_new_RETURN Pympf_new Pympf_new_PROTO;
static Pympf_dealloc_RETURN Pympf_dealloc Pympf_dealloc_PROTO;
static Pympf_convert_arg_RETURN Pympf_convert_arg Pympf_convert_arg_PROTO;

#if PY_MAJOR_VERSION < 3
#define export_gmpy(m) { \
      PyObject *d; \
      static void *Pygmpy_API[Pygmpy_API_pointers]; \
      PyObject *c_api_object; \
\
      Pygmpy_API[Pympz_Type_NUM] = (void*)&Pympz_Type;\
      Pygmpy_API[Pympq_Type_NUM] = (void*)&Pympq_Type;\
      Pygmpy_API[Pympf_Type_NUM] = (void*)&Pympf_Type;\
\
      Pygmpy_API[Pympz_new_NUM] = (void*)Pympz_new;\
      Pygmpy_API[Pympz_dealloc_NUM] = (void*)Pympz_dealloc;\
      Pygmpy_API[Pympz_convert_arg_NUM] = (void*)Pympz_convert_arg;\
      Pygmpy_API[Pympq_new_NUM] = (void*)Pympq_new;\
      Pygmpy_API[Pympq_dealloc_NUM] = (void*)Pympq_dealloc;\
      Pygmpy_API[Pympq_convert_arg_NUM] = (void*)Pympq_convert_arg;\
      Pygmpy_API[Pympf_new_NUM] = (void*)Pympf_new;\
      Pygmpy_API[Pympf_dealloc_NUM] = (void*)Pympf_dealloc;\
      Pygmpy_API[Pympf_convert_arg_NUM] = (void*)Pympf_convert_arg;\
\
      c_api_object = PyCObject_FromVoidPtr((void*)Pygmpy_API, NULL);\
      d = PyModule_GetDict(m);\
      PyDict_SetItemString(d, "_C_API", c_api_object);\
    }
#endif

#else
/* This section is used in other C-coded modules that use gmpy's API */

static void **Pygmpy_API;

#define Pympz_Check(op) \
  ((op)->ob_type == (PyTypeObject *)Pygmpy_API[Pympz_Type_NUM])
#define Pympz_Type (*(PyTypeObject *)Pygmpy_API[Pympz_Type_NUM])
#define Pympq_Check(op) \
  ((op)->ob_type == (PyTypeObject *)Pygmpy_API[Pympq_Type_NUM])
#define Pympq_Type (*(PyTypeObject *)Pygmpy_API[Pympq_Type_NUM])
#define Pympf_Check(op) \
  ((op)->ob_type == (PyTypeObject *)Pygmpy_API[Pympf_Type_NUM])
#define Pympf_Type (*(PyTypeObject *)Pygmpy_API[Pympf_Type_NUM])

#define Pympz_new \
  (*(Pympz_new_RETURN (*)Pympz_new_PROTO) Pygmpy_API[Pympz_new_NUM])
#define Pympz_dealloc \
  (*(Pympz_dealloc_RETURN (*)Pympz_dealloc_PROTO) Pygmpy_API[Pympz_dealloc_NUM])
#define Pympz_convert_arg \
  (*(Pympz_convert_arg_RETURN (*)Pympz_convert_arg_PROTO) Pygmpy_API[Pympz_convert_arg_NUM])
#define Pympq_new \
  (*(Pympq_new_RETURN (*)Pympq_new_PROTO) Pygmpy_API[Pympq_new_NUM])
#define Pympq_dealloc \
  (*(Pympq_dealloc_RETURN (*)Pympq_dealloc_PROTO) Pygmpy_API[Pympq_dealloc_NUM])
#define Pympq_convert_arg \
  (*(Pympq_convert_arg_RETURN (*)Pympq_convert_arg_PROTO) Pygmpy_API[Pympq_convert_arg_NUM])
#define Pympf_new \
  (*(Pympf_new_RETURN (*)Pympf_new_PROTO) Pygmpy_API[Pympf_new_NUM])
#define Pympf_dealloc \
  (*(Pympf_dealloc_RETURN (*)Pympf_dealloc_PROTO) Pygmpy_API[Pympf_dealloc_NUM])
#define Pympf_convert_arg \
  (*(Pympf_convert_arg_RETURN (*)Pympf_convert_arg_PROTO) Pygmpy_API[Pympf_convert_arg_NUM])

#if PY_MAJOR_VERSION < 3
#define import_gmpy() \
     {  \
       PyObject *module = PyImport_ImportModule("gmpy");\
       if (module != NULL) { \
         PyObject *module_dict = PyModule_GetDict(module); \
         PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API"); \
         if (PyCObject_Check(c_api_object)) { \
           Pygmpy_API = (void **)PyCObject_AsVoidPtr(c_api_object); \
         } \
       } \
     }
#else
static int
import_gmpy(void)
{
    Pygmpy_API = (void **)PyCapsule_Import("gmpy._C_API", 0);
    return (Pygmpy_API != NULL) ? 0 : -1;
}
#endif
#endif

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
