/* Low-level utility routines.
 *
 * The routines are:
 *   mpz_set_PyInt(mpz_t, PyObject)         Python 2.X only.
 *   mpz_set_PyLong(mpz_t, PyObject)
 *   mpf_normalize(mpf_t)
 *   Pympz_FROM_MPZ(mpz_t)
 *   Pympq_FROM_MPQ(mpq_t)
 *   Pympf_FROM_MPF(mpf_t)
 */

#ifdef FALSE
#if PY2
static void
mpz_set_PyInt(mpz_t rop, PyObject *obj)
{
    assert(PyInt_Check(obj));
    mpz_set_si(rop, PyInt_AsLong(obj));
    return;
}
#endif
#endif

static PympzObject *
Pympz_FROM_MPZ(mpz_t z)
{
    PympzObject *self;

    if(!(self = PyObject_New(PympzObject, &Pympz_Type)))
        return NULL;
    self->z[0] = z[0];
    return self;
}

#ifdef FALSE
static PympqObject *
Pympq_FROM_MPQ(mpq_t q)
{
    PympqObject *self;

    if(!(self = PyObject_New(PympqObject, &Pympq_Type)))
        return NULL;
    self->q[0] = q[0];
    return self;
}

static PympfObject *
Pympf_FROM_MPF(mpf_t f, unsigned int bits)
{
    PympfObject *self;

    if(!(self = PyObject_New(PympfObject, &Pympf_Type)))
        return NULL;
    if(bits < options.minprec)
        bits = options.minprec;
    self->f[0] = f[0];
    self->rebits = bits;
    return self;
}
#endif
/* End of low-level utility routines. */
