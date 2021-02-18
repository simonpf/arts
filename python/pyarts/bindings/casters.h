#define PY_SSIZE_T_CLEAN
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <mystring.h>
#include <iostream>

namespace py = pybind11;

namespace pybind11 { namespace detail {
    template <> struct type_caster<String> {
    public:
        PYBIND11_TYPE_CASTER(my_basic_string<char>, _("arts_string"));

        bool load(handle src, bool) {
            /* Extract PyObject from handle */
            std::string string =  src.cast<std::string>();
            value = string;
            return true;
        }

        /**
         * Conversion part 2 (C++ -> Python): convert an inty instance into
         * a Python object. The second and third arguments are used to
         * indicate the return value policy and parent object (for
         * ``return_value_policy::reference_internal``) and are generally
         * ignored by implicit casters.
         */
        static handle cast(String source,
                           return_value_policy /* policy */,
                           handle /* parent */) {
            return PyUnicode_FromString(source.c_str());
        }
    };
}}
