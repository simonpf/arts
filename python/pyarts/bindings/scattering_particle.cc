#define PY_SSIZE_T_CLEAN
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen_tensor.h>
#include <interactive_workspace.h>
#include <scattering/particle.h>
#include <casters.h>
#include <stdint.h>

namespace py = pybind11;

scattering::Particle&  from_wsv(intptr_t ptr) {
    return *reinterpret_cast<scattering::Particle *>(ptr);
}

void copy_to_wsv(intptr_t ws_ptr, Index id, const scattering::Particle &src) {
  auto workspace = reinterpret_cast<InteractiveWorkspace *>(ws_ptr);
  scattering::Particle &dst =
      *reinterpret_cast<scattering::Particle *>(workspace->operator[](id));
  dst = src;
}

void scattering_particle_module(py::module &bindings_module) {
  py::module m = bindings_module.def_submodule("scattering_particle");
  m.def("from_wsv", &from_wsv);
  m.def("copy_to_wsv", &copy_to_wsv);
  py::class_<scattering::Particle>(m, "Particle")
  // Constructors
  .def(py::init<>())
  .def(py::init<scattering::ParticleProperties,scattering::SingleScatteringData>())
  .def(py::init<double,double,double,scattering::SingleScatteringData>())
  .def(py::init<const scattering::Particle &>())
  // Class methods
  .def("operator=", &scattering::Particle::operator=)
  .def("copy", &scattering::Particle::copy)
  .def("get_name", &scattering::Particle::get_name)
  .def("get_source", &scattering::Particle::get_source)
  .def("get_refractive_index", &scattering::Particle::get_refractive_index)
  .def("get_particle_type", &scattering::Particle::get_particle_type)
  .def("get_data_format", &scattering::Particle::get_data_format)
  .def("get_mass", &scattering::Particle::get_mass)
  .def("get_d_max", &scattering::Particle::get_d_max)
  .def("get_d_eq", &scattering::Particle::get_d_eq)
  .def("get_d_aero", &scattering::Particle::get_d_aero)
  .def("get_f_grid", &scattering::Particle::get_f_grid)
  .def("get_t_grid", &scattering::Particle::get_t_grid)
  .def("get_lon_inc", &scattering::Particle::get_lon_inc)
  .def("get_lat_inc", &scattering::Particle::get_lat_inc)
  .def("get_lon_scat", &scattering::Particle::get_lon_scat)
  .def("get_lat_scat", &scattering::Particle::get_lat_scat)
  .def("interpolate_temperature", (scattering::SingleScatteringData(scattering::Particle::*)(double)const ) &scattering::Particle::interpolate_temperature)
  .def("to_spectral", &scattering::Particle::to_spectral)
  .def("to_lab_frame", (scattering::Particle(scattering::Particle::*)(Eigen::Index, Eigen::Index, Eigen::Index)const ) &scattering::Particle::to_lab_frame)
  .def("regrid", &scattering::Particle::regrid)
  .def("set_stokes_dim", &scattering::Particle::set_stokes_dim)
  .def("needs_t_interpolation", &scattering::Particle::needs_t_interpolation)
  .def("get_data", &scattering::Particle::get_data)
  .def("get_phase_function", &scattering::Particle::get_phase_function)
  .def("get_phase_function_spectral", &scattering::Particle::get_phase_function_spectral)
  .def("get_phase_matrix_data", &scattering::Particle::get_phase_matrix_data)
  .def("get_phase_matrix_data_spectral", &scattering::Particle::get_phase_matrix_data_spectral)
  .def("get_scattering_matrix", &scattering::Particle::get_scattering_matrix)
  .def("get_extinction_coeff", &scattering::Particle::get_extinction_coeff)
  .def("get_extinction_matrix_data", &scattering::Particle::get_extinction_matrix_data)
  .def("get_extinction_matrix", &scattering::Particle::get_extinction_matrix)
  .def("get_absorption_coeff", &scattering::Particle::get_absorption_coeff)
  .def("get_absorption_vector_data", &scattering::Particle::get_absorption_vector_data)
  .def("get_absorption_vector", &scattering::Particle::get_absorption_vector)
  .def("get_forward_scattering_coeff", &scattering::Particle::get_forward_scattering_coeff)
  .def("get_backward_scattering_coeff", &scattering::Particle::get_backward_scattering_coeff)
  // Data members
;
}
