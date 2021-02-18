#define PY_SSIZE_T_CLEAN
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <scattering/single_scattering_data.h>
#include <casters.h>

namespace py = pybind11;

void single_scattering_data_module(py::module &bindings_module) {

  py::module m = bindings_module.def_submodule("single_scattering_data");

  py::enum_<scattering::DataFormat>(m, "DataFormat")
      .value("Gridded", scattering::DataFormat::Gridded)
      .value("Spectral", scattering::DataFormat::Spectral)
      .value("FullySpectral", scattering::DataFormat::FullySpectral)
      .export_values();

  py::enum_<scattering::ParticleType>(m, "ParticleType")
      .value("Random", scattering::ParticleType::Random)
      .value("AzimuthallyRandom", scattering::ParticleType::AzimuthallyRandom)
      .value("General", scattering::ParticleType::General)
      .export_values();

  py::class_<scattering::SingleScatteringData>(m, "SingleScatteringData")
      // Constructors
      .def(py::init<>())
      .def(py::init<const scattering::SingleScatteringData &>())
      .def(py::init<scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Tensor<double, 7>,
                    scattering::eigen::Tensor<double, 7>,
                    scattering::eigen::Tensor<double, 7>,
                    scattering::eigen::Tensor<double, 7>,
                    scattering::eigen::Tensor<double, 7>>())
      .def(py::init<scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::ParticleType>())
      .def(py::init<scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::sht::SHT,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>>())
      .def(py::init<scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>,
                    scattering::eigen::Tensor<std::complex<double>, 6>>())
      .def(py::init<scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Vector<double>,
                    scattering::eigen::Index,
                    scattering::ParticleType>())
      // Class methods
      .def("copy", &scattering::SingleScatteringData::copy)
      .def("get_particle_type",
           &scattering::SingleScatteringData::get_particle_type)
      .def("get_data_format",
           &scattering::SingleScatteringData::get_data_format)
      .def("get_f_grid", &scattering::SingleScatteringData::get_f_grid)
      .def("get_t_grid", &scattering::SingleScatteringData::get_t_grid)
      .def("get_lon_inc", &scattering::SingleScatteringData::get_lon_inc)
      .def("get_lat_inc", &scattering::SingleScatteringData::get_lat_inc)
      .def("get_lon_scat", &scattering::SingleScatteringData::get_lon_scat)
      .def("get_lat_scat", &scattering::SingleScatteringData::get_lat_scat)
      .def("get_n_freqs", &scattering::SingleScatteringData::get_n_freqs)
      .def("get_n_temps", &scattering::SingleScatteringData::get_n_temps)
      .def("get_n_lon_inc", &scattering::SingleScatteringData::get_n_lon_inc)
      .def("get_n_lat_inc", &scattering::SingleScatteringData::get_n_lat_inc)
      .def("get_n_lon_scat", &scattering::SingleScatteringData::get_n_lon_scat)
      .def("get_n_lat_scat", &scattering::SingleScatteringData::get_n_lat_scat)
      .def("get_l_max_scat", &scattering::SingleScatteringData::get_l_max_scat)
      .def("get_m_max_scat", &scattering::SingleScatteringData::get_m_max_scat)
      .def("get_stokes_dim", &scattering::SingleScatteringData::get_stokes_dim)
      .def("set_data", &scattering::SingleScatteringData::set_data)
      .def("interpolate_frequency",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::eigen::Vector<double>) const) &
               scattering::SingleScatteringData::interpolate_frequency)
      .def("interpolate_temperature",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::eigen::Vector<double>, bool) const) &
               scattering::SingleScatteringData::interpolate_temperature)
      .def("interpolate_angles",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::eigen::Vector<double>,
               scattering::eigen::Vector<double>,
               scattering::eigen::Vector<double>,
               scattering::eigen::Vector<double>) const) &
               scattering::SingleScatteringData::interpolate_angles)
      .def("downsample_scattering_angles",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::eigen::Vector<double>,
               scattering::eigen::Vector<double>) const) &
               scattering::SingleScatteringData::downsample_scattering_angles)
      .def("get_phase_function",
           &scattering::SingleScatteringData::get_phase_function)
      .def("get_phase_function_spectral",
           &scattering::SingleScatteringData::get_phase_function_spectral)
      .def("get_scattering_matrix",
           &scattering::SingleScatteringData::get_scattering_matrix)
      .def("get_phase_matrix_data",
           &scattering::SingleScatteringData::get_phase_matrix_data)
      .def("get_phase_matrix_data_spectral",
           &scattering::SingleScatteringData::get_phase_matrix_data_spectral)
      .def("get_extinction_matrix_data",
           &scattering::SingleScatteringData::get_extinction_matrix_data)
      .def("get_extinction_coeff",
           &scattering::SingleScatteringData::get_extinction_coeff)
      .def("get_extinction_matrix",
           &scattering::SingleScatteringData::get_extinction_matrix)
      .def("get_absorption_vector_data",
           &scattering::SingleScatteringData::get_absorption_vector_data)
      .def("get_absorption_coeff",
           &scattering::SingleScatteringData::get_absorption_coeff)
      .def("get_absorption_vector",
           &scattering::SingleScatteringData::get_absorption_vector)
      .def("get_forward_scattering_coeff",
           &scattering::SingleScatteringData::get_forward_scattering_coeff)
      .def("get_backward_scattering_coeff",
           &scattering::SingleScatteringData::get_backward_scattering_coeff)
      .def("__iadd__", &scattering::SingleScatteringData::operator+=)
      .def("__add__", &scattering::SingleScatteringData::operator+)
      .def("__imul__", &scattering::SingleScatteringData::operator*=)
      .def("__mul__", &scattering::SingleScatteringData::operator*)
      .def("normalize", &scattering::SingleScatteringData::normalize)
      .def("regrid", &scattering::SingleScatteringData::regrid)
      .def("set_number_of_scattering_coeffs",
           &scattering::SingleScatteringData::set_number_of_scattering_coeffs)
      .def("to_gridded", &scattering::SingleScatteringData::to_gridded)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)() const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::eigen::Index, scattering::eigen::Index) const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(scattering::eigen::Index,
                                                    scattering::eigen::Index,
                                                    scattering::eigen::Index,
                                                    scattering::eigen::Index)
                const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_lab_frame",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(scattering::eigen::Index,
                                                    scattering::eigen::Index,
                                                    scattering::eigen::Index)
                const) &
               scattering::SingleScatteringData::to_lab_frame)
      .def("set_stokes_dim", &scattering::SingleScatteringData::set_stokes_dim)
      // Data members
      ;
}
