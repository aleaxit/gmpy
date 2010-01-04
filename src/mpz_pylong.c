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
#if PyLong_SHIFT >= GMP_NUMB_BITS
#error "Python limb larger than GMP limb !!!"
#endif

#ifndef ABS
#define ABS(a)  (((a) < 0) ? -(a) : (a))
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
size_t
mpn_sizebits(mp_ptr up, size_t un) {
  size_t cnt;
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
size_t
pylong_sizebits(digit *digits, size_t size) {
  size_t cnt;
  digit x;
  if (size==0) return 0;
  cnt = (size - 1) * PyLong_SHIFT;
  x = digits[size - 1];
#if PyLong_SHIFT > 32
  if ((x >> 32) != 0) { x >>= 32; cnt += 32; }
#endif
#if PyLong_SHIFT > 16
  if ((x >> 16) != 0) { x >>= 16; cnt += 16; }
#endif
#if PyLong_SHIFT > 8
  if ((x >>  8) != 0) { x >>=  8; cnt += 8; }
#endif
  return cnt + ((x & 0x80) ? 8 : __sizebits_tab[x]);
}

/* mpn -> pylong conversion */

size_t
mpn_pylong_size (mp_ptr up, size_t un)
{
  return (mpn_sizebits(up, un) + PyLong_SHIFT - 1) / PyLong_SHIFT;
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
  ssize_t bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (un == 0) {
    while (size) digits[--size]=0;
    return;
  }

  i = un - 1;
  n1 = up[i];
  bit_pos = size * PyLong_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= PyLong_SHIFT;
      while (bit_pos >= 0)
        {
          *--s = (n1 >> bit_pos) & PyLong_MASK;
          bit_pos -= PyLong_SHIFT;
        }
      if (i == 0)
        break;
      n0 = (n1 << -bit_pos) & PyLong_MASK;
      n1 = up[--i];
      bit_pos += GMP_NUMB_BITS;
      *--s = (digit)(n0 | (n1 >> bit_pos));
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
  size_t i;
  ssize_t bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (size == 0) {
    while (un) up[--un]=0;
    return;
  }

  i = un - 1;
  n1 = 0;
  bit_pos = size * PyLong_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= PyLong_SHIFT;
      while (bit_pos >= 0)
        {
          d = (mp_limb_t) *--s;
          n1 |= (d << bit_pos) & GMP_NUMB_MASK;
          bit_pos -= PyLong_SHIFT;
        }
      if (i == 0)
        break;
      d = (mp_limb_t) *--s;
      /* add some high bits of d; maybe none if bit_pos=-SHIFT */
      up[i--] = n1 | (d & PyLong_MASK) >> -bit_pos;
      bit_pos += GMP_NUMB_BITS;
      n1 = (d << bit_pos) & GMP_NUMB_MASK;
    }
  up[0] = n1;
}


/************************************************************/

/* Hashing functions */

#define LONG_BIT_SHIFT  (8*sizeof(long) - PyLong_SHIFT)

/*
 * for an mpz, this number has to be multiplied by the sign
 * also remember to catch -1 and map it to -2 !
 */
long
mpn_pythonhash (mp_ptr up, mp_size_t un)
{
  mp_limb_t n1, n0;
  mp_size_t i;
  ssize_t bit_pos;
  long x = 0;

  /* input length 0 is special ! */
  if (un == 0) return 0;

  i = un - 1;
  n1 = up[i];
  {
    size_t bits;
    bits = mpn_sizebits(up, un) + PyLong_SHIFT - 1;
    bits -= bits % PyLong_SHIFT;
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
          x = ((x << PyLong_SHIFT) & ~(long)PyLong_MASK) | ((x >> LONG_BIT_SHIFT) & (long)PyLong_MASK);
          /* Shifting to the right by more than wordsize bits
             actually shifts by (wordsize % 32) bits -- which is
             *not* the intended behavior here. */
      if (bit_pos <= 8*sizeof(mp_limb_t))
            x += (n1 >> bit_pos) & (long)PyLong_MASK;
          bit_pos -= PyLong_SHIFT;
        }
      i--;
      if (i < 0)
        break;
      n0 = (n1 << -bit_pos) & (long)PyLong_MASK;
      n1 = up[i];
      bit_pos += GMP_NUMB_BITS;
      /* Force a native long #-bits (32 or 64) circular shift */
      x = ((x << PyLong_SHIFT) & ~(long)PyLong_MASK) | ((x >> LONG_BIT_SHIFT) & (long)PyLong_MASK);
      x += (long)(n0 | (n1 >> bit_pos));
      bit_pos -= PyLong_SHIFT;
    }

  return x;
}


/* mpz python hash */
long
mpz_pythonhash(mpz_srcptr z)
{
  long x = mpn_pythonhash(z->_mp_d, ABS(z->_mp_size));
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
  size_t size = mpn_pylong_size(z->_mp_d, ABS(z->_mp_size));
  PyLongObject *lptr = PyObject_NEW_VAR(PyLongObject, &PyLong_Type, size);

  if (lptr != NULL)
  {
    mpn_get_pylong(lptr->ob_digit, size, z->_mp_d, ABS(z->_mp_size));
    if (z->_mp_size < 0)
      Py_SIZE(lptr) = -(Py_SIZE(lptr));
  }

  return (PyObject *) lptr;
}

/* pylong -> mpz conversion */
int
mpz_set_PyLong(mpz_ptr z, PyObject * lsrc)
{
  register PyLongObject * lptr = (PyLongObject *) lsrc;
  ssize_t size;

  if (lptr==NULL || !PyLong_Check(lptr)) {
    PyErr_BadInternalCall();
    return -1;
  }

  size = (ssize_t)mpn_size_from_pylong(lptr->ob_digit, ABS(Py_SIZE(lptr)));

  if (z->_mp_alloc < size)
    _mpz_realloc (z, (mp_size_t)size);

  mpn_set_pylong(z->_mp_d, size, lptr->ob_digit, ABS(Py_SIZE(lptr)));
  z->_mp_size = (int)(Py_SIZE(lptr) < 0 ? -size : size);

  return (int)size;
}

