/* mpz <-> pylong conversion and "pythonhash" for mpz
 *
 * Originally written for sage (http://sagemath.org) by Gonzalo Tornari­a
 * <tornaria@math.utexas.edu>. If you improve on these functions, please
 * contribute them back to sage by posting to sage-devel@googlegroups.com
 * or by sending an email to the original author.
 *
 * Integration with gmpy by Case Van Horsen <casevh@gmail.com>.
 *
 * License: LGPL v2 or later
 *
 */

/* This file created by merging mpn_pylong and mpz_pylong. Permission
 * was granted by the original author to make this code available under
 * the LGPLv2+ license.
 */

/* This code assumes that SHIFT < GMP_NUMB_BITS */
#if Py2or3Long_SHIFT >= GMP_NUMB_BITS
#error "Python limb larger than GMP limb !!!"
#endif

/* Use these "portable" (I hope) sizebits functions.
 * We could implement this in terms of count_leading_zeros from GMP,
 * but it is not exported!
 */
static const
unsigned char
__sizebits_tab[128] =
{
  0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};

#if GMP_LIMB_BITS > 64
#error "word size > 64 unsupported"
#endif

static inline
unsigned long
mpn_sizebits(mp_ptr up, mp_size_t un) {
  unsigned long cnt;
  mp_limb_t x;
  if (un==0) return 0;
  cnt = (un - 1) * GMP_NUMB_BITS;
  x = up[un - 1];
#if GMP_LIMB_BITS > 32
  if ((x >> 32) != 0) { x >>= 32; cnt += 32; }
#endif
#if GMP_LIMB_BITS > 16
  if ((x >> 16) != 0) { x >>= 16; cnt += 16; }
#endif
#if GMP_LIMB_BITS > 8
  if ((x >>  8) != 0) { x >>=  8; cnt += 8; }
#endif
  return cnt + ((x & 0x80) ? 8 : __sizebits_tab[x]);
}

static inline
unsigned long
pylong_sizebits(digit *digits, size_t size) {
  unsigned long cnt;
  digit x;
  if (size==0) return 0;
  cnt = (size - 1) * Py2or3Long_SHIFT;
  x = digits[size - 1];
#if Py2or3Long_SHIFT > 32
  if ((x >> 32) != 0) { x >>= 32; cnt += 32; }
#endif
#if Py2or3Long_SHIFT > 16
  if ((x >> 16) != 0) { x >>= 16; cnt += 16; }
#endif
#if Py2or3Long_SHIFT > 8
  if ((x >>  8) != 0) { x >>=  8; cnt += 8; }
#endif
  return cnt + ((x & 0x80) ? 8 : __sizebits_tab[x]);
}

/* mpn -> pylong conversion */

int
mpn_pylong_size (mp_ptr up, size_t un)
{
  return (mpn_sizebits(up, un) + Py2or3Long_SHIFT - 1) / Py2or3Long_SHIFT;
}

/* this is based from GMP code in mpn/get_str.c */

/* Assume digits points to a chunk of size size
 * where size >= mpn_pylong_size(up, un)
 */
void
mpn_get_pylong (digit *digits, size_t size, mp_ptr up, size_t un)
{
  mp_limb_t n1, n0;
  size_t i;
  int bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (un == 0) {
    while (size) digits[--size]=0;
    return;
  }

  i = un - 1;
  n1 = up[i];
  bit_pos = size * Py2or3Long_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= Py2or3Long_SHIFT;
      while (bit_pos >= 0)
        {
          *--s = (n1 >> bit_pos) & Py2or3Long_MASK;
          bit_pos -= Py2or3Long_SHIFT;
        }
      if (i == 0)
        break;
      n0 = (n1 << -bit_pos) & Py2or3Long_MASK;
      n1 = up[--i];
      bit_pos += GMP_NUMB_BITS;
      *--s = n0 | (n1 >> bit_pos);
    }
}

/* pylong -> mpn conversion */

size_t
mpn_size_from_pylong (digit *digits, size_t size)
{
  return (pylong_sizebits(digits, size) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
}

void
mpn_set_pylong (mp_ptr up, size_t un, digit *digits, size_t size)
{
  mp_limb_t n1, d;
  mp_size_t i;
  int bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (size == 0) {
    while (un) up[--un]=0;
    return;
  }

  i = un - 1;
  n1 = 0;
  bit_pos = size * Py2or3Long_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= Py2or3Long_SHIFT;
      while (bit_pos >= 0)
        {
          d = (mp_limb_t) *--s;
          n1 |= (d << bit_pos) & GMP_NUMB_MASK;
          bit_pos -= Py2or3Long_SHIFT;
        }
      if (i == 0)
        break;
      d = (mp_limb_t) *--s;
      /* add some high bits of d; maybe none if bit_pos=-SHIFT */
      up[i--] = n1 | (d & Py2or3Long_MASK) >> -bit_pos;
      bit_pos += GMP_NUMB_BITS;
      n1 = (d << bit_pos) & GMP_NUMB_MASK;
    }
  up[0] = n1;
}


/************************************************************/

/* Hashing functions */

#define LONG_BIT_SHIFT  (8*sizeof(long) - Py2or3Long_SHIFT)

/*
 * for an mpz, this number has to be multiplied by the sign
 * also remember to catch -1 and map it to -2 !
 */
long
mpn_pythonhash (mp_ptr up, mp_size_t un)
{
  mp_limb_t n1, n0;
  mp_size_t i;
  int bit_pos;
  long x = 0;

  /* input length 0 is special ! */
  if (un == 0) return 0;

  i = un - 1;
  n1 = up[i];
  {
    unsigned long bits;
    bits = mpn_sizebits(up, un) + Py2or3Long_SHIFT - 1;
    bits -= bits % Py2or3Long_SHIFT;
    /* position of the MSW in base 2^SHIFT, counted from the MSW in
     * the GMP representation (in base 2^GMP_NUMB_BITS)
     */
    bit_pos = bits - i * GMP_NUMB_BITS;
  }

  for (;;)
    {
      while (bit_pos >= 0)
        {
          /* Force a native long #-bits (32 or 64) circular shift */
          x = ((x << Py2or3Long_SHIFT) & ~(long)Py2or3Long_MASK) | ((x >> LONG_BIT_SHIFT) & (long)Py2or3Long_MASK);
	  /* Shifting to the right by more than wordsize bits
             actually shifts by (wordsize % 32) bits -- which is
             *not* the intended behavior here. */
	  if (bit_pos <= 8*sizeof(mp_limb_t))
            x += (n1 >> bit_pos) & (long)Py2or3Long_MASK;
          bit_pos -= Py2or3Long_SHIFT;
        }
      i--;
      if (i < 0)
        break;
      n0 = (n1 << -bit_pos) & (long)Py2or3Long_MASK;
      n1 = up[i];
      bit_pos += GMP_NUMB_BITS;
      /* Force a native long #-bits (32 or 64) circular shift */
      x = ((x << Py2or3Long_SHIFT) & ~(long)Py2or3Long_MASK) | ((x >> LONG_BIT_SHIFT) & (long)Py2or3Long_MASK);
      x += n0 | (n1 >> bit_pos);
      bit_pos -= Py2or3Long_SHIFT;
    }

  return x;
}


/* mpz python hash */
long
mpz_pythonhash(mpz_srcptr z)
{
  long x = mpn_pythonhash(z->_mp_d, abs(z->_mp_size));
  if (z->_mp_size < 0)
    x = -x;
  if (x == -1)
    x = -2;
  return x;
}

/* mpz -> pylong conversion */
PyObject *
mpz_get_PyLong(mpz_srcptr z)
{
  size_t size = mpn_pylong_size(z->_mp_d, abs(z->_mp_size));
  PyLongObject *lptr = PyObject_NEW_VAR(PyLongObject, &PyLong_Type, size);

  if (lptr != NULL)
  {
    mpn_get_pylong(lptr->ob_digit, size, z->_mp_d, abs(z->_mp_size));
    if (z->_mp_size < 0)
#if PY_MAJOR_VERSION >= 3
      lptr->ob_base.ob_size = -(lptr->ob_base.ob_size);
#else
      lptr->ob_size = -(lptr->ob_size);
#endif
  }

  return (PyObject *) lptr;
}

/* pylong -> mpz conversion */
int
mpz_set_PyLong(mpz_ptr z, PyObject * lsrc)
{
  register PyLongObject * lptr = (PyLongObject *) lsrc;
  size_t size;
  //~ int i;

  if (lptr==NULL || !PyLong_Check(lptr)) {
    PyErr_BadInternalCall();
    return -1;
  }

#if PY_MAJOR_VERSION >= 3
  size = mpn_size_from_pylong(lptr->ob_digit, abs(lptr->ob_base.ob_size));
#else
  size = mpn_size_from_pylong(lptr->ob_digit, abs(lptr->ob_size));
#endif

  if (z->_mp_alloc < size)
    _mpz_realloc (z, size);

#if PY_MAJOR_VERSION >= 3
  mpn_set_pylong(z->_mp_d, size, lptr->ob_digit, abs(lptr->ob_base.ob_size));
  z->_mp_size = lptr->ob_base.ob_size < 0 ? -size : size;
#else
  mpn_set_pylong(z->_mp_d, size, lptr->ob_digit, abs(lptr->ob_size));
  z->_mp_size = lptr->ob_size < 0 ? -size : size;
#endif

  return size;
}

