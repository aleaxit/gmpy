/* utility macros for creating the correct return type */


/*
 * Create two 'mpz' and a 2-tuple.
 */

#define CREATE_TWO_MPZ_TUPLE(q, r, t)\
    q = Pympz_new();\
    r = Pympz_new();\
    t = PyTuple_New(2);\
    if(!q || !r || !t) {\
        Py_XDECREF(t);\
        Py_XDECREF((PyObject*)q);\
        Py_XDECREF((PyObject*)r);\
        return NULL;\
    }

/* utility macros for argument parsing */

/*
 * Verify that a function has only one argument, and convert that argument
 * to a C-long. Only applies to gmpy.fname()."msg" should be an error
 * message that includes the function name. Replaces ONE_ARG.
 */

#define PARSE_ONE_CLONG(var, msg)\
    if (PyTuple_GET_SIZE(args) != 1) {\
        PyErr_SetString(PyExc_TypeError, msg);\
        return NULL;\
    } else {\
        *var = clong_From_Integer(PyTuple_GET_ITEM(args, 0)); \
        if(*var == -1 && PyErr_Occurred()) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
    }

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpz. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces SELF_MPZ_NO_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ(msg) \
    if(self && CHECK_MPZANY(self)) {\
        if (PyTuple_GET_SIZE(args) != 0) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        Py_INCREF(self);\
    } else {\
        if (PyTuple_GET_SIZE(args) != 1) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        self = PyTuple_GET_ITEM(args, 0);\
        if(CHECK_MPZANY(self)) {\
            Py_INCREF((PyObject*)self);\
        } else {\
            self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
        }\
        if(!self) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
    }

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpf. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces SELF_MPZ_NO_ARG.
 */

#define PARSE_ONE_MPF(msg) \
    if(self && Pympf_Check(self)) {\
        if (PyTuple_GET_SIZE(args) != 0) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        Py_INCREF(self);\
    } else {\
        if (PyTuple_GET_SIZE(args) != 1) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        self = PyTuple_GET_ITEM(args, 0);\
        if(Pympf_Check(self)) {\
            Py_INCREF((PyObject*)self);\
        } else {\
            self = (PyObject*)Pympf_From_Float(PyTuple_GET_ITEM(args, 0), 0);\
        }\
        if(!self) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
    }

/*
 * Parses one argument into "self" and an optional second argument into
 * 'var". The second argument is converted into a C long. If there is not a
 * second argument, "var" is unchanged. Is faster, but not as generic, as
 * using PyArg_ParseTuple with "|l". It supports either gmpy.fname(z,l) or
 * z.fname(l). "self" must be decref'ed. "var" must be a pointer to a long.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces some uses of SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_OPT_CLONG(var, msg) \
    if(self && CHECK_MPZANY(self)) {\
        if (PyTuple_GET_SIZE(args) == 1) {\
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if(*var == -1 && PyErr_Occurred()) {\
                PyErr_SetString(PyExc_TypeError, msg);\
                return NULL;\
            }\
        } else if (PyTuple_GET_SIZE(args) > 1) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        Py_INCREF(self);\
    } else {\
        if (PyTuple_GET_SIZE(args) == 2) {\
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if(*var == -1 && PyErr_Occurred()) {\
                PyErr_SetString(PyExc_TypeError, msg);\
                return NULL;\
            }\
            self = PyTuple_GET_ITEM(args, 0);\
            if(CHECK_MPZANY(self)) {\
                Py_INCREF((PyObject*)self);\
            } else {\
                self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
            }\
        } else if (PyTuple_GET_SIZE(args) == 1) {\
            self = PyTuple_GET_ITEM(args, 0);\
            if(CHECK_MPZANY(self)) {\
                Py_INCREF((PyObject*)self);\
            } else {\
                self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
            }\
        } else {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        if(!self) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
    }

/*
 * Parses one argument into "self" and a required second argument into
 * 'var". The second argument is converted into a C long. Is faster, but not
 * as generic, as using PyArg_ParseTuple with "l". It supports either
 * gmpy.fname(z,l) or z.fname(l). "self" must be decref'ed. "var" must be a
 * pointer to a long. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces some uses of
 * SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_REQ_CLONG(var, msg) \
    if(self && CHECK_MPZANY(self)) {\
        if (PyTuple_GET_SIZE(args) != 1) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        } else {\
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if(*var == -1 && PyErr_Occurred()) {\
                PyErr_SetString(PyExc_TypeError, msg);\
                return NULL;\
            }\
        }\
        Py_INCREF(self);\
    } else {\
        if (PyTuple_GET_SIZE(args) != 2) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        } else {\
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if(*var == -1 && PyErr_Occurred()) {\
                PyErr_SetString(PyExc_TypeError, msg);\
                return NULL;\
            }\
            self = PyTuple_GET_ITEM(args, 0);\
            if(CHECK_MPZANY(self)) {\
                Py_INCREF((PyObject*)self);\
            } else {\
                self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
            }\
        }\
        if(!self) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
    }

/*
 * Parses two, and only two, arguments into "self" and "var" and converts
 * them both to mpz. Is faster, but not as generic, as using PyArg_ParseTuple.
 * It supports either gmpy.fname(z,z) or z.fname(z). "self" & "var" must be
 * decref'ed after use. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces
 * SELF_MPZ_ONE_ARG_CONVERTED(var).
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_TWO_MPZ(var, msg) \
    if(self && CHECK_MPZANY(self)) {\
        if (PyTuple_GET_SIZE(args) != 1) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        var = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
        if(!var) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        Py_INCREF(self);\
    } else {\
        if (PyTuple_GET_SIZE(args) != 2) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            return NULL;\
        }\
        self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));\
        var = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));\
        if(!self || !var) {\
            PyErr_SetString(PyExc_TypeError, msg);\
            Py_XDECREF((PyObject*)self);\
            Py_XDECREF((PyObject*)var);\
            return NULL;\
        }\
    }

/* Define three different versions of the SELF_NO_ARG macro. Under Python
   2.x, self is NULL when a function is called via gmpy.fname(..). But
   under Python 3.x, self is a module. */

#define SELF_MPQ_NO_ARG \
    if(self && Pympq_Check(self)) { \
        if(!PyArg_ParseTuple(args, "")) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&", Pympq_convert_arg, &self)) \
            return NULL; \
    }
#define SELF_MPF_NO_ARG \
    if(self && Pympf_Check(self)) { \
        if(!PyArg_ParseTuple(args, "")) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&", Pympf_convert_arg, &self)) \
            return NULL; \
    }

#define SELF_MPQ_ONE_ARG(fm, var) \
    if(self && Pympq_Check(self)) { \
        if(!PyArg_ParseTuple(args, fm, var)) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&" fm, Pympq_convert_arg, &self, var)) \
            return NULL; \
    }

#define SELF_MPF_ONE_ARG(fm, var) \
    if(self && Pympf_Check(self)) { \
        if(!PyArg_ParseTuple(args, fm, var)) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&" fm, Pympf_convert_arg, &self, var)) \
            return NULL; \
    }

#define SELF_MPQ_ONE_ARG_CONVERTED(var) \
    if(self && Pympq_Check(self)) { \
        if(args && !PyArg_ParseTuple(args, "O&", Pympq_convert_arg, var)) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&O&", Pympq_convert_arg,&self, \
                Pympq_convert_arg,var)) \
            return NULL; \
    }
#define SELF_MPF_ONE_ARG_CONVERTED(var) \
    if(self && Pympf_Check(self)) { \
        if(args && !PyArg_ParseTuple(args, "O&", Pympf_convert_arg, var)) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&O&", Pympf_convert_arg,&self, \
                Pympf_convert_arg,var)) \
            return NULL; \
    }


#define SELF_MPF_ONE_ARG_CONVERTED_OPT(var) \
    if(self && Pympf_Check(self)) { \
        if(args && !PyArg_ParseTuple(args, "|O&", Pympf_convert_arg,var)) \
            return NULL; \
        Py_INCREF(self); \
    } else { \
        if(!PyArg_ParseTuple(args, "O&|O&", Pympf_convert_arg,&self, \
                Pympf_convert_arg,var)) \
            return NULL; \
    }


#define TWO_ARG_CONVERTED(converter, var1, var2) \
    if(!PyArg_ParseTuple(args, "O&O&", converter,var1, converter,var2)) \
        return NULL;

