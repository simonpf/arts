#define PY_SSIZE_T_CLEAN
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <interactive_workspace.h>
#include <iostream>
#include <casters.h>

namespace py = pybind11;


void single_scattering_data_module(py::module &);
void scattering_particle_module(py::module &);
void xml_io_module(py::module &bindings_module);
void messages_module(py::module &bindings_module);

void print(const String& s) {
    std::cout << s << std::endl;
}

PYBIND11_MODULE(cxx, m) {
    single_scattering_data_module(m);
    scattering_particle_module(m);
    xml_io_module(m);
    messages_module(m);
    m.def("print", &print);
}
