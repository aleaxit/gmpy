/* gmpy_mpf.c
 *
 * Functions that operate strictly on mpf.
 *
 * This file should be considered part of gmpy2.c
 */

static PyObject *
Pympf_getprec2(PympfObject *self, void *closure)
{
    return PyIntOrLong_FromSize_t((size_t)mpfr_get_prec(self->f));
}

