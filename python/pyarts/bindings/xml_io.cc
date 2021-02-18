#define PY_SSIZE_T_CLEAN
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <scattering/single_scattering_data.h>
#include <casters.h>

#include "xml_io.h"
#include "xml_io_types.h"

namespace py = pybind11;

void xml_io_module(py::module &bindings_module) {
    py::module m = bindings_module.def_submodule("xml_io");

    py::enum_<FileType>(m, "FileType")
        .value("FILE_TYPE_ASCII", FileType::FILE_TYPE_ASCII)
        .value("FILE_TYPE_ZIPPED_ASCII", FileType::FILE_TYPE_ZIPPED_ASCII)
        .value("FILE_TYPE_BINARY", FileType::FILE_TYPE_BINARY)
        .export_values();

    m.def("xml_read_from_file", (void (*)(const String&,
                                          ScatteringParticle&,
                                          const Verbosity&)) &xml_read_from_file);
    m.def("xml_write_to_file", (void (*)(const String&,
                                         const ScatteringParticle &,
                                         const FileType,
                                         const Index,
                                         const Verbosity&)) &xml_write_to_file);
}
