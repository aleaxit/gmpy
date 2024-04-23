/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_demo.c                                                            *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_demo contains a collection of demonstration programs based on     *
 * C code originally provided in the /demos/ directory of the GMP          *
 * source distribution. gmpy2_demo utilizes the gmpy2's C-API to interact  *
 * with the gmpy2.                                                         *
 *                                                                         *
 * Since the original C code in the /demos/ directory is licensed under    *
 * the GPL (in contrast to to GMP which is licensed under LGPL), this file *
 * is also licensed under GPL.                                             *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003 Pearu Peterson                         *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017 Case Van Horsen                              *
 *                                                                         *
 * GMPY2_DEMO is free software: you can redistribute it and/or modify it   *
 * under the terms of the GNU General Public License as published by the   *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2_DEMO is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * General Public License for more details.                                *
 *                                                                         *
 * You should have received a copy of the GNU General Public License along *
 * with GMPY2_DEMO; if not, see <http://www.gnu.org/licenses/>             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* gmpy2_demo.c is an updated version of pysymbolicext.c originally written
 * by Pearu Peterson. The original license is below.
 */

/* PySymbolic GMP extensions (in connection with gmpy).

   1) Factoring with Pollard's rho method.
   Copied relevant functions from demos/factorice.c of the GMP distribution,
   and modified them to be used in Python.

Copyright 2000 Pearu Peterson all rights reserved,
Pearu Peterson <pearu@ioc.ee>
Permission to use, modify, and distribute this software is given under the
terms of the LGPL.  See http://www.fsf.org

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
$Revision: 1.2 $
$Date: 2003/08/08 08:57:05 $
Pearu Peterson

Bug fix for more than 32 iterations in pollard_rho method.
  (Patch courtesy rel...@osagesoftware.com)

*/
#include "Python.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gmpy2.h>

#define MPZ(obj) (((MPZ_Object*)(obj))->z)

static unsigned add[] = {4, 2, 4, 2, 4, 6, 2, 6};

#if defined (__hpux) || defined (__alpha)  || defined (__svr4__) || defined (__SVR4)
/* HPUX lacks random().  DEC OSF/1 1.2 random() returns a double.  */
long mrand48 ();
static long
random ()
{
  return mrand48 ();
}
#else
/* Glibc stdlib.h has "int32_t random();" which, on i386 at least, conflicts
   with a redeclaration as "long". */
#if defined(_MSC_VER)
long random() { return rand(); }
#endif
#ifndef __GLIBC__
long random ();
#endif
#endif

static void
res_append_si(PyObject *res,signed long int f,unsigned long int c) {
  if (c) {
    PyObject *pair = (PyObject *)PyTuple_New(2);
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set_si(MPZ(z),f);
      PyTuple_SetItem(pair,0,z);
    }
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set_ui(MPZ(z),c);
      PyTuple_SetItem(pair,1,z);
    }
    PyList_Append(res,pair);
    Py_DECREF(pair);
  }
}

static void
res_append(PyObject *res,unsigned long int f,unsigned long int c) {
  if (c) {
    PyObject *pair = (PyObject *)PyTuple_New(2);
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set_ui(MPZ(z),f);
      PyTuple_SetItem(pair,0,z);
    }
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set_ui(MPZ(z),c);
      PyTuple_SetItem(pair,1,z);
    }
    PyList_Append(res,pair);
    Py_DECREF(pair);
  }
}

static void
res_append_mpz(PyObject *res,mpz_t f,unsigned long int c) {
  if (c) {
    PyObject *pair = (PyObject *)PyTuple_New(2);
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set(MPZ(z),f);
      PyTuple_SetItem(pair,0,z);
    }
    {
      PyObject *z = (PyObject *)GMPy_MPZ_New(NULL);
      mpz_set_ui(MPZ(z),c);
      PyTuple_SetItem(pair,1,z);
    }
    PyList_Append(res,pair);
    Py_DECREF(pair);
  }
}

static void
factor_using_division (mpz_t t, unsigned int limit, PyObject *res)
{
  mpz_t q, r;
  unsigned long int f;
  int ai;
  unsigned *addv = add;
  unsigned int failures;

  unsigned long int count;

  mpz_init (q);
  mpz_init (r);


  count = 0;
  f = mpz_scan1 (t, 0);
  mpz_fdiv_q_2exp (t, t, f);
  res_append(res,2,f);

  count = 0;
  for (;;)
    {
      mpz_tdiv_qr_ui (q, r, t, 3);
      if (mpz_cmp_ui (r, 0) != 0)
        break;
      mpz_set (t, q);
      count++;
    }
  res_append(res,3,count);
  count = 0;
  for (;;)
    {
      mpz_tdiv_qr_ui (q, r, t, 5);
      if (mpz_cmp_ui (r, 0) != 0)
        break;
      mpz_set (t, q);
      count++;
    }
  res_append(res,5,count);

  failures = 0;
  f = 7;
  ai = 0;
  count = 0;
  while (mpz_cmp_ui (t, 1) != 0)
    {
      mpz_tdiv_qr_ui (q, r, t, f);
      if (mpz_cmp_ui (r, 0) != 0)
        {
          res_append(res,f,count);
          count = 0;
          f += addv[ai];
          if (mpz_cmp_ui (q, f) < 0)
            break;
          ai = (ai + 1) & 7;
          failures++;
          if (failures > limit)
            break;
        }
      else
        {
          mpz_swap (t, q);
          failures = 0;
          count++;
        }
    }
  res_append(res,f,count);
  mpz_clear (q);
  mpz_clear (r);
}

static void
factor_using_division_2kp (mpz_t t, unsigned int limit, unsigned long p, PyObject *res)
{
  mpz_t r;
  mpz_t f;
  unsigned int k;
  unsigned long int count;
  mpz_init (r);
  mpz_init_set_ui (f, 2 * p);
  mpz_add_ui (f, f, 1);
  for (k = 1; k < limit; k++)
    {
      mpz_tdiv_r (r, t, f);
      count = 0;
      while (mpz_cmp_ui (r, 0) == 0)
        {
          mpz_tdiv_q (t, t, f);
          mpz_tdiv_r (r, t, f);
          count++;
        }
      res_append_mpz(res,f,count);
      mpz_add_ui (f, f, 2 * p);
    }

  mpz_clear (f);
  mpz_clear (r);
}

static
void
factor_using_pollard_rho (mpz_t n, int a_int, unsigned long p,PyObject *res)
{
  mpz_t x, x1, y, P;
  mpz_t a;
  mpz_t g;
  mpz_t t1, t2;
  mpz_t kz, lz, iz;
  int c;
  unsigned long int count;

  mpz_init (g);
  mpz_init (t1);
  mpz_init (t2);

  mpz_init_set_si (a, a_int);
  mpz_init_set_si (y, 2);
  mpz_init_set_si (x, 2);
  mpz_init_set_si (x1, 2);

  mpz_init (iz);
  mpz_init_set_si (kz, 1);
  mpz_init_set_si (lz, 1);
  mpz_init_set_ui (P, 1);
  c = 0;
  count = 0;
  while (mpz_cmp_ui (n, 1) != 0)
    {
S2:
      if (p != 0)
        {
          mpz_powm_ui (x, x, p, n); mpz_add (x, x, a);
        }
      else
        {
          mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);
        }
      mpz_sub (t1, x1, x); mpz_mul (t2, P, t1); mpz_mod (P, t2, n);
      c++;
      if (c == 20)
        {
          c = 0;
          mpz_gcd (g, P, n);
          if (mpz_cmp_ui (g, 1) != 0)
            goto S4;
          mpz_set (y, x);
        }
/* S3: */
      mpz_sub_ui (kz, kz, 1);
      if (mpz_cmp_ui(kz,0) > 0)
        goto S2;

      mpz_gcd (g, P, n);
      if (mpz_cmp_ui (g, 1) != 0)
        goto S4;

      mpz_set (x1, x);
      mpz_set (kz, lz);
      mpz_mul_ui (lz, lz, 2);

      // for loop with integer index works fine for k < 2**31
      // using mpz_t allows unlimited range
      for ( mpz_set (iz, kz);  mpz_cmp_ui(iz,0) > 0; mpz_sub_ui (iz, iz, 1) )
        {
          if (p != 0)
            {
              mpz_powm_ui (x, x, p, n); mpz_add (x, x, a);
            }
          else
            {
              mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);
            }
        }
      mpz_set (y, x);
      c = 0;
      goto S2;
S4:
      do
        {
          if (p != 0)
            {
              mpz_powm_ui (y, y, p, n); mpz_add (y, y, a);
            }
          else
            {
              mpz_mul (y, y, y); mpz_add (y, y, a); mpz_mod (y, y, n);
            }
          mpz_sub (t1, x1, y); mpz_gcd (g, t1, n);
        }
      while (mpz_cmp_ui (g, 1) == 0);

      if (!mpz_probab_prime_p (g, 3))
        {
          do
            a_int = random ();
          while (a_int == -2 || a_int == 0);
          factor_using_pollard_rho (g, a_int, p, res);
          break;
        }
      else
        count ++;

      res_append_mpz(res,g,count);
      count = 0;
      mpz_fdiv_q (n, n, g);
      mpz_mod (x, x, n);
      mpz_mod (x1, x1, n);
      mpz_mod (y, y, n);
      if (mpz_probab_prime_p (n, 3))
        {
          count++;
          break;
        }
    }
  res_append_mpz(res,n,count);
  mpz_clear (iz);
  mpz_clear (kz);
  mpz_clear (lz);
  mpz_clear (g);
  mpz_clear (P);
  mpz_clear (t2);
  mpz_clear (t1);
  mpz_clear (a);
  mpz_clear (x1);
  mpz_clear (x);
  mpz_clear (y);
}

static void
factor (mpz_t t, unsigned long p,PyObject *res)
{
  unsigned int division_limit;

  /* Set the trial division limit according the size of t.  */
  division_limit = mpz_sizeinbase (t, 2);
  if (division_limit > 1000)
    division_limit = 1000 * 1000;
  else
    division_limit = division_limit * division_limit;

  if (p != 0)
    factor_using_division_2kp (t, division_limit / 10, p, res);
  else
    factor_using_division (t, division_limit, res);

  if (mpz_cmp_ui (t, 1) != 0) {
    if (mpz_probab_prime_p (t, 3))
      res_append_mpz(res,t,1);
    else
      factor_using_pollard_rho (t, 1, p, res);
  }
}

static char doc_factor[] =
"factor(t,m=0) -> prime factors of t (modulo m)\n\
\n\
Prime decomposition of t (modulo m).\n\
factor(t,m) returns a list of tuples (f,p) where p is the number\n\
of prime factors f in term t. p is always positive.\n\
t can be also zero or negative.\n\
For m=0 the following condition holds\n\
     t == reduce(lambda r,pm:r*pm[0]**pm[1],factor(t),1L).";
static PyObject *
Pygmpy2_demo_factor(PyObject *self, PyObject *args)
{
  mpz_t t;
  unsigned long p;
  PyObject *t_py = NULL;
  PyObject *p_py = NULL;
  PyObject *res = NULL;
  if (!PyArg_ParseTuple(args, "O&|O&",\
                        GMPy_MPZ_ConvertArg,&t_py,\
                        GMPy_MPZ_ConvertArg,&p_py))
    return NULL;

  res = PyList_New(0);
  if (p_py==NULL)
    p = 0;
  else
    p = mpz_get_ui(MPZ(p_py));
  mpz_init_set(t,MPZ(t_py));

  if (mpz_sgn(t)==0) {
    res_append(res,0,1);
    return res;
  }
  if (mpz_sgn(t)==-1) {
    res_append_si(res,-1,1);
    mpz_neg(t,t);
  }

  factor(t,p,res);

  if (PyList_Size(res)==0)
    res_append(res,1,1);
  return res;
}

static PyMethodDef Pygmpy2_demo_methods [] =
{
    { "factor", Pygmpy2_demo_factor, METH_VARARGS, doc_factor },
    { NULL, NULL}
};

#define INITERROR return NULL
static struct PyModuleDef Pygmpy2_demo_module = {
    PyModuleDef_HEAD_INIT,
    "gmpy2_demo",
    NULL,
    -1,
    Pygmpy2_demo_methods,
    NULL,
    NULL,
    NULL,
    NULL
};
#ifdef _MSC_VER
__declspec(dllexport)
#endif
PyObject *
PyInit_gmpy2_demo(void)
{
    PyObject *gmpy2_demo_module = NULL;
    gmpy2_demo_module = PyModule_Create(&Pygmpy2_demo_module);

    import_gmpy2();

    return gmpy2_demo_module;
}
