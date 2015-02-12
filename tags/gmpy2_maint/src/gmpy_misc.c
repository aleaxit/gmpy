/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_misc.c                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015 Case Van Horsen                                          *
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

/* Miscellaneous module-level functions and helper functions. */

PyDoc_STRVAR(doc_license,
"license() -> string\n\n"
"Return string giving license information.");

static PyObject *
Pygmpy_get_license(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_license);
}

PyDoc_STRVAR(doc_version,
"version() -> string\n\n"
"Return string giving current GMPY2 version.");

static PyObject *
Pygmpy_get_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_version);
}

PyDoc_STRVAR(doc_cvsid,
"_cvsid() -> string\n\n"
"Return string giving current GMPY2 cvs Id.");

static PyObject *
Pygmpy_get_cvsid(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", _gmpy_cvs);
}

PyDoc_STRVAR(doc_mp_version,
"mp_version() -> string\n\n"
"Return string giving the name and version of the multiple precision\n"
"library used.");

static PyObject *
Pygmpy_get_mp_version(PyObject *self, PyObject *args)
{
#ifndef __MPIR_VERSION
    return Py2or3String_FromFormat("%s %s", "GMP", gmp_version);
#else
    return Py2or3String_FromFormat("%s %s", "MPIR", mpir_version);
#endif
}

PyDoc_STRVAR(doc_mpfr_version,
"mpfr_version() -> string\n\n"
"Return string giving current MPFR version. Return None if MPFR\n"
"support is not available.");

static PyObject *
Pygmpy_get_mpfr_version(PyObject *self, PyObject *args)
{
#ifdef WITHMPFR
    return Py2or3String_FromFormat("%s %s", "MPFR",
                                   MPFR_VERSION_STRING);
#else
    Py_RETURN_NONE;
#endif
}

PyDoc_STRVAR(doc_mpc_version,
"mpc_version() -> string\n\n"
"Return string giving current MPC version. Return None if MPC\n"
"support is not available.");

static PyObject *
Pygmpy_get_mpc_version(PyObject *self, PyObject *args)
{
#ifdef WITHMPC
    return Py2or3String_FromFormat("%s %s", "MPC",
                                   MPC_VERSION_STRING);
#else
    Py_RETURN_NONE;
#endif
}

PyDoc_STRVAR(doc_mp_limbsize,
"mp_limbsize() -> integer\n\n\
Return the number of bits per limb.");

static PyObject *
Pygmpy_get_mp_limbsize(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", mp_bits_per_limb);
}

/*
 * access cache options
 */

PyDoc_STRVAR(doc_get_cache,
"get_cache() -> (cache_size, object_size)\n\n\
Return the current cache size (number of objects) and maximum size\n\
per object (number of limbs) for all GMPY2 objects.");

static PyObject *
Pygmpy_get_cache(PyObject *self, PyObject *args)
{
    return Py_BuildValue("(ii)", global.cache_size, global.cache_obsize);
}

PyDoc_STRVAR(doc_set_cache,
"set_cache(cache_size, object_size)\n\n\
Set the current cache size (number of objects) and the maximum size\n\
per object (number of limbs). Raises ValueError if cache size exceeds\n\
1000 or object size exceeds 16384.");

static PyObject *
Pygmpy_set_cache(PyObject *self, PyObject *args)
{
    int newcache = -1, newsize = -1;

    if (!PyArg_ParseTuple(args, "ii", &newcache, &newsize))
        return NULL;
    if (newcache<0 || newcache>MAX_CACHE) {
        VALUE_ERROR("cache size must between 0 and 1000");
        return NULL;
    }
    if (newsize<0 || newsize>MAX_CACHE_LIMBS) {
        VALUE_ERROR("object size must between 0 and 16384");
        return NULL;
    }

    global.cache_size = newcache;
    global.cache_obsize = newsize;
    set_zcache();
    set_pympzcache();
    set_pympqcache();
    set_pyxmpzcache();
#ifdef WITHMPFR
    set_pympfrcache();
#endif
#ifdef WITHMPC
    set_pympccache();
#endif
    Py_RETURN_NONE;
}
