/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_prp.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2011 David Cleaver                                            *
 *                                                                         *
 * Copyright 2012 Case Van Horsen                                          *
 *                                                                         *
 * The original file is available at:                                      *
 *   <http://sourceforge.net/projects/mpzprp/files/>                       *
 *                                                                         *
 * Modified by Case Van Horsen for inclusion into GMPY2.                   *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ******************************************************************
 * mpz_prp: (also called a Fermat pseudoprime)
 * A "pseudoprime" to the base a is a composite number n such that,
 * (a,n)=1 and a^(n-1) = 1 mod n
 * ******************************************************************/
 
PyDoc_STRVAR(doc_mpz_is_fermat_prp,
"is_fermat_prp(n,a) -> boolean\n\n"
"Return True if 'n' is a Fermat pseudoprime to the base 'a'.");

static PyObject *
GMPY_mpz_is_fermat_prp(PyObject *self, PyObject *args)
{
    PympzObject *a, *n;
    PyObject *result = 0;
    mpz_t res, nm1;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_fermat_prp() requires 2 integer arguments");
        return NULL;
    }
    
    n = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    a = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (!a || !n) {
        TYPE_ERROR("is_fermat_prp() requires 2 integer arguments");
        goto cleanup;
    }
    
    /* Take advantage of the cache of mpz_t objects maintained by GMPY2 to
     * avoid memory allocations. */
     
    mpz_inoc(res);
    mpz_inoc(nm1);
           
    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_fermat_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }
    
    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_fermat_prp() requires 'n' be greater than 0");
        goto cleanup;
    }
    
    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,n) */
    mpz_gcd(res, n->z, a->z);
    if (mpz_cmp_ui(res, 1) > 0) {
        VALUE_ERROR("is_fermat_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }

    mpz_set(nm1, n->z);
    mpz_sub_ui(nm1, nm1, 1);
    mpz_powm(res, a->z, nm1, n->z);

    if (mpz_cmp_ui(res, 1) == 0)
        result = Py_True;
    else
        result = Py_False;
        
  cleanup:
    Py_XINCREF(result);
    mpz_cloc(res);
    mpz_cloc(nm1);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *************************************************************************
 * mpz_euler_prp: (also called a Solovay-Strassen pseudoprime)
 * An "Euler pseudoprime" to the base a is an odd composite number n with,
 * (a,n)=1 such that a^((n-1)/2)=(a/n) mod n [(a/n) is the Jacobi symbol]
 * *************************************************************************/
 
PyDoc_STRVAR(doc_mpz_is_euler_prp,
"is_euler_prp(n,a) -> boolean\n\n"
"Return True if 'n' is an Euler (also known as Solovay-Strassen)\n"
"pseudoprime to the base 'a'.");

static PyObject *
GMPY_mpz_is_euler_prp(PyObject *self, PyObject *args)
{
    PympzObject *a, *n;
    PyObject *result = 0;
    mpz_t res, exp;
    int ret;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_euler_prp() requires 2 integer arguments");
        return NULL;
    }
    
    n = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    a = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (!a || !n) {
        TYPE_ERROR("is_euler_prp() requires 2 integer arguments");
        goto cleanup;
    }
    
    /* Take advantage of the cache of mpz_t objects maintained by GMPY2 to
     * avoid memory allocations. */
     
    mpz_inoc(res);
    mpz_inoc(exp);
           
    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_euler_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }
    
    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_euler_prp() requires 'n' be greater than 0");
        goto cleanup;
    }
    
    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,b) */
    mpz_gcd(res, n->z, a->z);
    if (mpz_cmp_ui(res, 1) > 0) {
        VALUE_ERROR("is_euler_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }

    mpz_set(exp, n->z);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(res, a->z, exp, n->z);

    /* reuse exp to calculate jacobi(a,n) mod n */
    ret = mpz_jacobi(a->z,n->z);
    mpz_set(exp, n->z);
    if (ret == -1)
        mpz_sub_ui(exp, exp, 1);
    else if (ret == 1)
        mpz_add_ui(exp, exp, 1);
    mpz_mod(exp, exp, n->z);

    if (mpz_cmp(res, exp) == 0)
        result = Py_True;
    else
        result = Py_False;
        
  cleanup:
    Py_XINCREF(result);
    mpz_cloc(res);
    mpz_cloc(exp);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *********************************************************************************************
 * mpz_sprp: (also called a Miller-Rabin pseudoprime)
 * A "strong pseudoprime" to the base a is an odd composite n = (2^r)*s+1 with s odd such that
 * either a^s == 1 mod n, or a^((2^t)*s) == -1 mod n, for some integer t, with 0 <= t < r.
 * *********************************************************************************************/
 
PyDoc_STRVAR(doc_mpz_is_strong_prp,
"is_strong_prp(n,a) -> boolean\n\n"
"Return True if 'n' is an strong (also known as Miller-Rabin)\n"
"pseudoprime to the base 'a'.");

static PyObject *
GMPY_mpz_is_strong_prp(PyObject *self, PyObject *args)
{
    PympzObject *a, *n;
    PyObject *result = 0;
    mpz_t s, nm1, mpz_test;
    mpir_ui r = 0;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_strong_prp() requires 2 integer arguments");
        return NULL;
    }
    
    n = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    a = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (!a || !n) {
        TYPE_ERROR("is_strong_prp() requires 2 integer arguments");
        goto cleanup;
    }
    
    /* Take advantage of the cache of mpz_t objects maintained by GMPY2 to
     * avoid memory allocations. */
     
    mpz_inoc(s);
    mpz_inoc(nm1);
    mpz_inoc(mpz_test);
           
    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_strong_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }
    
    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_strong_prp() requires 'n' be greater than 0");
        goto cleanup;
    }
    
    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,b) */
    mpz_gcd(s, n->z, a->z);
    if (mpz_cmp_ui(s, 1) > 0) {
        VALUE_ERROR("is_strong_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }
    
    mpz_set(nm1, n->z);
    mpz_sub_ui(nm1, nm1, 1);

    /* Find s and r satisfying: n-1=(2^r)*s, s odd */
    r = mpz_scan1(nm1, 0);
    mpz_fdiv_q_2exp(s, nm1, r);


    /* Check a^((2^t)*s) mod n for 0 <= t < r */
    mpz_powm(mpz_test, a->z, s, n->z);
    if ((mpz_cmp_ui(mpz_test, 1) == 0) || (mpz_cmp(mpz_test, nm1) == 0)) {
        result = Py_True;
        goto cleanup;
    }

    while (--r) {
        /* mpz_test = mpz_test^2%n */
        mpz_mul(mpz_test, mpz_test, mpz_test);
        mpz_mod(mpz_test, mpz_test, n->z);

        if (mpz_cmp(mpz_test, nm1) == 0) {
            result = Py_True;
            goto cleanup;
        }
    }

    result = Py_False;
  cleanup:
    Py_XINCREF(result);
    mpz_cloc(s);
    mpz_cloc(nm1);
    mpz_cloc(mpz_test);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *************************************************************************
 * mpz_fibonacci_prp:
 * A "Fibonacci pseudoprime" with parameters (P,Q), P > 0, Q=+/-1, is a
 * composite n for which V_n == P mod n
 * [V is the Lucas V sequence with parameters P,Q]
 * *************************************************************************/
 
PyDoc_STRVAR(doc_mpz_is_fibonacci_prp,
"is_fibonacci_prp(n,p,q) -> boolean\n\n"
"Return True if 'n' is an Fibonacci pseudoprime with parameters (p,q).\n"
"p > 0; q = +/-1; True if lucasv(p,q,n) = p (mod n).");

static PyObject *
GMPY_mpz_is_fibonacci_prp(PyObject *self, PyObject *args)
{
    PympzObject *n, *p, *q;
    PyObject *result = 0;
    mpz_t pmodn, zP;
    /* used for calculating the Lucas V sequence */
    mpz_t vl, vh, ql, qh, tmp;
    size_t s = 0, j = 0;
    
    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("is_fibonacci_prp() requires 3 integer arguments");
        return NULL;
    }
    
    /* Take advantage of the cache of mpz_t objects maintained by GMPY2 to
     * avoid memory allocations. */
    
    mpz_inoc(pmodn);
    mpz_inoc(zP);
    mpz_inoc(vl);        
    mpz_inoc(vh);
    mpz_inoc(ql);
    mpz_inoc(qh);
    mpz_inoc(tmp);
    
    n = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    p = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
    p = Pympz_From_Integer(PyTuple_GET_ITEM(args, 2));
    if (!n || !p || !q) {
        TYPE_ERROR("is_fibonacci_prp() requires 3 integer arguments");
        goto cleanup;
    }
    
    /* Check if p*p - 4*q == 0. */
    
    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in is_fibonacci_prp()");
        goto cleanup;
    }

    /* Verify q = +/-1 */
    
    if ((mpz_cmp_si(q->z, 1) && mpz_cmp_si(q->z, -1)) || (mpz_sgn(p->z) <= 0)) {
        VALUE_ERROR("invalid values for p,q in is_fibonacci_prp()");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_fibonacci_prp() requires 'n' be greater than 0");
        goto cleanup;
    }
    
    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    mpz_set_ui(zP, p->z);
    mpz_mod(pmodn, zP, n->z);

    /* mpz_lucasvmod(res, p, q, n, n); */
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    s = mpz_scan1(n->z, 0);
    for (j = mpz_sizeinbase(n->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n);
        if (mpz_tstbit(n->z,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul(qh, ql, q->z);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p->z);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    for (j = 1; j <= s; j++) {
        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n->z);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n->z);
    }
    
    /* vl contains our return value */
    mpz_mod(vl, vl, n->z); 

    if (mpz_cmp(vl, pmodn) == 0) 
        result = Py_True;
    else
        result = Py_False;

  cleanup:
    Py_XINCREF(result);
    mpz_cloc(pmodn);
    mpz_cloc(zP);
    mpz_cloc(vl);        
    mpz_cloc(vh);
    mpz_cloc(ql);
    mpz_cloc(qh);
    mpz_cloc(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)n);
    return result;
}


#if 0
/* *******************************************************************************
 * mpz_lucas_prp:
 * A "Lucas pseudoprime" with parameters (P,Q) is a composite n with D=P^2-4Q,
 * (n,2QD)=1 such that U_(n-(D/n)) == 0 mod n [(D/n) is the Jacobi symbol]
 * *******************************************************************************/
int mpz_lucas_prp(mpz_t n, long int p, long int q)
{
  mpz_t zD;
  mpz_t res;
  mpz_t index;
  mpz_t uh, vl, vh, ql, qh, tmp; /* used for calculating the Lucas U sequence */
  int s = 0, j = 0;
  int ret = 0;
  long int d = p*p - 4*q;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init(index);
  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_si(res, zD, q);
  mpz_mul_ui(res, res, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    mpz_clear(index);
    return PRP_COMPOSITE;
  }

  /* index = n-(D/n), where (D/n) is the Jacobi symbol */
  mpz_set(index, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(index, index, 1);
  else if (ret == 1)
    mpz_sub_ui(index, index, 1);

  /* mpz_lucasumod(res, p, q, index, n); */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  s = mpz_scan1(index, 0);
  for (j = mpz_sizeinbase(index,2)-1; j >= s+1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(index,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  for (j = 1; j <= s; j++)
  {
    /* uh = uh*vl (mod n) */
    mpz_mul(uh, uh, vl);
    mpz_mod(uh, uh, n);

    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);
  }

  mpz_mod(res, uh, n); /* uh contains our return value */

  mpz_clear(zD);
  mpz_clear(index);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);

  if (mpz_cmp_ui(res, 0) == 0)
  {
    mpz_clear(res);
    return PRP_PRP;
  }
  else
  {
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

}/* method mpz_lucas_prp */


/* *********************************************************************************************
 * mpz_stronglucas_prp:
 * A "strong Lucas pseudoprime" with parameters (P,Q) is a composite n = (2^r)*s+(D/n), where
 * s is odd, D=P^2-4Q, and (n,2QD)=1 such that either U_s == 0 mod n or V_((2^t)*s) == 0 mod n
 * for some t, 0 <= t < r. [(D/n) is the Jacobi symbol]
 * *********************************************************************************************/
int mpz_stronglucas_prp(mpz_t n, long int p, long int q)
{
  mpz_t zD;
  mpz_t s;
  mpz_t nmj; /* n minus jacobi(D/n) */
  mpz_t res;
  mpz_t uh, vl, vh, ql, qh, tmp; /* these are needed for the LucasU and LucasV part of this function */
  long int d = p*p - 4*q;
  unsigned long int r = 0;
  int ret = 0;
  int j = 0;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_si(res, zD, q);
  mpz_mul_ui(res, res, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init(s);
  mpz_init(nmj);

  /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
  mpz_set(nmj, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(nmj, nmj, 1);
  else if (ret == 1)
    mpz_sub_ui(nmj, nmj, 1);

  r = mpz_scan1(nmj, 0);
  mpz_fdiv_q_2exp(s, nmj, r);

  /* make sure U_s == 0 mod n or V_((2^t)*s) == 0 mod n, for some t, 0 <= t < r */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(s,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  mpz_mod(uh, uh, n);
  mpz_mod(vl, vl, n);

  /* uh contains LucasU_s and vl contains LucasV_s */
  if ((mpz_cmp_ui(uh, 0) == 0) || (mpz_cmp_ui(vl, 0) == 0))
  {
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    return PRP_PRP;
  }

  for (j = 1; j < r; j++)
  {
    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);

    if (mpz_cmp_ui(vl, 0) == 0)
    {
      mpz_clear(zD);
      mpz_clear(s);
      mpz_clear(nmj);
      mpz_clear(res);
      mpz_clear(uh);
      mpz_clear(vl);
      mpz_clear(vh);
      mpz_clear(ql);
      mpz_clear(qh);
      mpz_clear(tmp);
      return PRP_PRP;
    }
  }

  mpz_clear(zD);
  mpz_clear(s);
  mpz_clear(nmj);
  mpz_clear(res);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);
  return PRP_COMPOSITE;

}/* method mpz_stronglucas_prp */


/* *******************************************************************************************
 * mpz_extrastronglucas_prp:
 * Let U_n = LucasU(p,1), V_n = LucasV(p,1), and D=p^2-4.
 * An "extra strong Lucas pseudoprime" to the base p is a composite n = (2^r)*s+(D/n), where
 * s is odd and (n,2D)=1, such that either U_s == 0 mod n or V_s == +/-2 mod n, or
 * V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1 [(D/n) is the Jacobi symbol]
 * *******************************************************************************************/
int mpz_extrastronglucas_prp(mpz_t n, long int p)
{
  mpz_t zD;
  mpz_t s;
  mpz_t nmj; /* n minus jacobi(D/n) */
  mpz_t res;
  mpz_t uh, vl, vh, ql, qh, tmp; /* these are needed for the LucasU and LucasV part of this function */
  long int d = p*p - 4;
  long int q = 1;
  unsigned long int r = 0;
  int ret = 0;
  int j = 0;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_ui(res, zD, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init(s);
  mpz_init(nmj);

  /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
  mpz_set(nmj, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(nmj, nmj, 1);
  else if (ret == 1)
    mpz_sub_ui(nmj, nmj, 1);

  r = mpz_scan1(nmj, 0);
  mpz_fdiv_q_2exp(s, nmj, r);

  /* make sure that either U_s == 0 mod n or V_s == +/-2 mod n, or */
  /* V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1           */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(s,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  mpz_mod(uh, uh, n);
  mpz_mod(vl, vl, n);

  /* uh contains LucasU_s and vl contains LucasV_s */
  if (   (mpz_cmp_ui(uh, 0) == 0) || (mpz_cmp_ui(vl, 0) == 0)
      || (mpz_cmp_si(vl, -2) == 0) || (mpz_cmp_si(vl, 2) == 0))
  {
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    return PRP_PRP;
  }

  for (j = 1; j < r-1; j++)
  {
    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);

    if (mpz_cmp_ui(vl, 0) == 0)
    {
      mpz_clear(zD);
      mpz_clear(s);
      mpz_clear(nmj);
      mpz_clear(res);
      mpz_clear(uh);
      mpz_clear(vl);
      mpz_clear(vh);
      mpz_clear(ql);
      mpz_clear(qh);
      mpz_clear(tmp);
      return PRP_PRP;
    }
  }

  mpz_clear(zD);
  mpz_clear(s);
  mpz_clear(nmj);
  mpz_clear(res);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);
  return PRP_COMPOSITE;

}/* method mpz_extrastronglucas_prp */


/* ***********************************************************************************************
 * mpz_selfridge_prp:
 * A "Lucas-Selfridge pseudoprime" n is a "Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the Lucas pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * ***********************************************************************************************/
int mpz_selfridge_prp(mpz_t n)
{
  long int d = 5, p = 1, q = 0;
  int max_d = 1000000;
  int jacobi = 0;
  mpz_t zD;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(zD, d);

  while (1)
  {
    jacobi = mpz_jacobi(zD, n);

    /* if jacobi == 0, d is a factor of n, therefore n is composite... */
    /* if d == n, then either n is either prime or 9... */
    if (jacobi == 0)
    {
      if ((mpz_cmpabs(zD, n) == 0) && (mpz_cmp_ui(zD, 9) != 0))
      {
        mpz_clear(zD);
        return PRP_PRIME;
      }
      else
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }
    if (jacobi == -1)
      break;

    /* if we get to the 5th d, make sure we aren't dealing with a square... */
    if (d == 13)
    {
      if (mpz_perfect_square_p(n))
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }

    if (d < 0)
    {
      d *= -1;
      d += 2;
    }
    else
    {
      d += 2;
      d *= -1;
    }

    /* make sure we don't search forever */
    if (d >= max_d)
    {
      mpz_clear(zD);
      return PRP_ERROR;
    }

    mpz_set_si(zD, d);
  }
  mpz_clear(zD);

  q = (1-d)/4;

  return mpz_lucas_prp(n, p, q);

}/* method mpz_selfridge_prp */


/* *********************************************************************************************************
 * mpz_strongselfridge_prp:
 * A "strong Lucas-Selfridge pseudoprime" n is a "strong Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the strong Lucase pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * **********************************************************************************************************/
int mpz_strongselfridge_prp(mpz_t n)
{
  long int d = 5, p = 1, q = 0;
  int max_d = 1000000;
  int jacobi = 0;
  mpz_t zD;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(zD, d);

  while (1)
  {
    jacobi = mpz_jacobi(zD, n);

    /* if jacobi == 0, d is a factor of n, therefore n is composite... */
    /* if d == n, then either n is either prime or 9... */
    if (jacobi == 0)
    {
      if ((mpz_cmpabs(zD, n) == 0) && (mpz_cmp_ui(zD, 9) != 0))
      {
        mpz_clear(zD);
        return PRP_PRIME;
      }
      else
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }
    if (jacobi == -1)
      break;

    /* if we get to the 5th d, make sure we aren't dealing with a square... */
    if (d == 13)
    {
      if (mpz_perfect_square_p(n))
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }

    if (d < 0)
    {
      d *= -1;
      d += 2;
    }
    else
    {
      d += 2;
      d *= -1;
    }

    /* make sure we don't search forever */
    if (d >= max_d)
    {
      mpz_clear(zD);
      return PRP_ERROR;
    }

    mpz_set_si(zD, d);
  }
  mpz_clear(zD);

  q = (1-d)/4;

  return mpz_stronglucas_prp(n, p, q);

}/* method mpz_strongselfridge_prp */


/* **********************************************************************************
 * mpz_bpsw_prp:
 * A "Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a Lucas pseudoprime using the Selfridge parameters.
 * **********************************************************************************/
int mpz_bpsw_prp(mpz_t n)
{
  int ret = 0;
  mpz_t two;

  mpz_init_set_ui(two, 2);

  ret = mpz_sprp(n, two);
  mpz_clear(two);

  /* with a base of 2, mpz_sprp, won't return PRP_ERROR */
  /* so, only check for PRP_COMPOSITE or PRP_PRIME here */
  if ((ret == PRP_COMPOSITE) || (ret == PRP_PRIME))
    return ret;

  return mpz_selfridge_prp(n);

}/* method mpz_bpsw_prp */


/* ****************************************************************************************
 * mpz_strongbpsw_prp:
 * A "strong Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a strong Lucas pseudoprime using the Selfridge parameters.
 * ****************************************************************************************/
int mpz_strongbpsw_prp(mpz_t n)
{
  int ret = 0;
  mpz_t two;

  mpz_init_set_ui(two, 2);

  ret = mpz_sprp(n, two);
  mpz_clear(two);

  /* with a base of 2, mpz_sprp, won't return PRP_ERROR */
  /* so, only check for PRP_COMPOSITE or PRP_PRIME here */
  if ((ret == PRP_COMPOSITE) || (ret == PRP_PRIME))
    return ret;

  return mpz_strongselfridge_prp(n);

}/* method mpz_strongbpsw_prp */

#endif
