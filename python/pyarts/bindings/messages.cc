#define PY_SSIZE_T_CLEAN
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <casters.h>

#include "messages.h"

namespace py = pybind11;

void messages_module(py::module &bindings_module) {
  py::module m = bindings_module.def_submodule("messages");
  py::class_<Verbosity>(m, "Verbosity")
      // Constructors
      .def(py::init<>())
      .def(py::init<Index, Index, Index>())
      // Class methods
      .def("valid", &Verbosity::valid)
      .def("get_agenda_verbosity", &Verbosity::get_agenda_verbosity)
      .def("get_screen_verbosity", &Verbosity::get_screen_verbosity)
      .def("get_file_verbosity", &Verbosity::get_file_verbosity)
      .def("is_main_agenda", &Verbosity::is_main_agenda)
      .def("set_agenda_verbosity", &Verbosity::set_agenda_verbosity)
      .def("set_screen_verbosity", &Verbosity::set_screen_verbosity)
      .def("set_file_verbosity", &Verbosity::set_file_verbosity)
      .def("set_main_agenda", &Verbosity::set_main_agenda)
      // Data members
      ;
}
