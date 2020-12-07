/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, rite to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   scattering_habit.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-15

  \brief  Implementation of scattering_habit.h
*/
#include "scattering_habit.h"
#include "auto_md.h"
#include <iomanip>

const String SCATSPECIES_MAINTAG = "Scattering species";

namespace detail {

/** Extract backward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scattering format.
 * @return Eigen tensor containing the backscattering coefficient in
 *     scattering-compatible format.
 */
EigenTensor<7> extract_backward_scattering_coeff(
    const EigenTensor<7> &phase_matrix) {
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(dimensions[5] - 1).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}

/** Extract forward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scattering format.
 * @return Eigen tensor containing the forward scattering coefficient in
 *     scattering-compatible format.
 */
EigenTensor<7> extract_forward_scattering_coeff(
    const EigenTensor<7> &phase_matrix) {
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(0).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}

/** Convert ARTS to scattering format.
 *
 * This function converts ARTS phase matrix data to scattering format.
 *
 * @param tensor ARTS phase matrix data.
 * @return Eigen tensor containing the phase matrix in scattering-compatible format.
 */
EigenTensor<7> arts_to_scattering(const Tensor7 &tensor) {
    EigenTensor<7> tensor_eigen = to_eigen(tensor);
    std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 5, 4, 3, 2, 6};
    EigenTensor<7> result = tensor_eigen.shuffle(shuffle_dimensions);
    return tensor_eigen.shuffle(shuffle_dimensions);
}

/** Convert ARTS to scattering format.
 *
 * This function converts a ARTS extinction matrix or absorption
 * vector data to scattering format
 *
 * @param tensor The input data to convert.
 * @return Eigen tensor containing the data in scattering-compatible format.
 */
EigenTensor<7> arts_to_scattering(const Tensor5 &tensor) {
  EigenTensor<5> tensor_eigen = to_eigen(tensor);
  auto result = scattering::eigen::unsqueeze<4, 5>(tensor_eigen);
  std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 3, 2, 4, 5, 6};
  return result.shuffle(shuffle_dimensions);
}

/** Convert ARTS legacy scattering data to new scattering format.
 *
 * @param arts_data SingleScatteringData object containing the scattering
 * data in legacy format.
 */
scattering::SingleScatteringData arts_to_scattering(
    const SingleScatteringData &arts_data) {
  EigenVector f_grid = to_eigen(arts_data.f_grid);
  EigenVector t_grid = to_eigen(arts_data.T_grid);
  EigenVector angle_grid_inc = EigenVector::Constant(1, M_PI);
  angle_grid_inc[0] = 0.0;
  EigenVector lon_scat = to_eigen(arts_data.aa_grid);
  if (lon_scat.size() == 0) {
      lon_scat = angle_grid_inc;
  }
  EigenVector lat_scat = to_eigen(arts_data.za_grid);
  lat_scat *= (M_PI / 180.0);
  auto phase_matrix = arts_to_scattering(arts_data.pha_mat_data);
  auto extinction_matrix = arts_to_scattering(arts_data.ext_mat_data);
  auto absorption_vector = arts_to_scattering(arts_data.abs_vec_data);
  auto backward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);
  auto forward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);

  auto dimensions = phase_matrix.dimensions();

  assert(phase_matrix.dimension(0) == f_grid.size());
  assert(phase_matrix.dimension(1) == t_grid.size());
  assert(phase_matrix.dimension(4) == lon_scat.size());
  assert(phase_matrix.dimension(5) == lat_scat.size());

  return scattering::SingleScatteringData(f_grid,
                                       t_grid,
                                       angle_grid_inc,
                                       angle_grid_inc,
                                       lon_scat,
                                       lat_scat,
                                       phase_matrix,
                                       extinction_matrix,
                                       absorption_vector,
                                       backward_scattering_coeff,
                                       forward_scattering_coeff);
}

}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
// ScatteringParticle and ArrayOfScatteringParticle
////////////////////////////////////////////////////////////////////////////////

std::ostream &operator<<(std::ostream &output,
                         const scattering::DataFormat format) {
  switch (format) {
    case scattering::DataFormat::Gridded:
      output << "gridded";
      break;
    case scattering::DataFormat::Spectral:
      output << "spectral";
      break;
  }
  return output;
}

std::ostream &operator<<(std::ostream &output,
                         const scattering::ParticleType format) {
  switch (format) {
    case scattering::ParticleType::Random:
      output << "random";
      break;
    case scattering::ParticleType::AzimuthallyRandom:
      output << "azimuthally random";
      break;
    case scattering::ParticleType::General:
      output << "general";
      break;
  }
  return output;
}

std::ostream &operator<<(std::ostream &output,
                         const ScatteringParticle &particle) {
  auto name = particle.get_name();
  auto source = particle.get_source();
  auto refractive_index = particle.get_refractive_index();

  output << "Scattering particle: " << name;
  if (source != "") {
    output << " (" << source << ")";
  }
  output << std::endl;
  if (refractive_index != "") {
    output << std::setw(30)
           << "\tRefractive index: " << particle.get_refractive_index()
           << std::endl;
  }
  output << std::setw(30) << "\tParticle type: " << particle.get_particle_type() << std::endl;
  output << std::setw(30) << "\tScattering data format: " << particle.get_data_format()
         << std::endl;
  output << std::setw(30) << "\tVolume-equivalent diameter: " << particle.get_d_eq() << std::endl;
  output << std::setw(30) << "\tMaximum diameter: " << particle.get_d_max() << std::endl;
  output << std::setw(30) << "\tMass: " << particle.get_d_max() << std::endl;
  output << std::endl;
  return output;
}

std::ostream &operator<<(std::ostream &output,
                         const ArrayOfScatteringParticle &array) {

    auto n_particles = array.nelem();
    auto p0 = array[0];
    auto pn = array[n_particles - 1];
    output << "Array of Scattering particles with " << n_particles;
    output << " particles:" << std::endl;
    output << std::setw(30) << "\tVolume-equivalent diameters: " << p0.get_d_eq()
           << ", ..., " << pn.get_d_eq() << std::endl;
    output << std::setw(30) << "\tMasses: " << p0.get_d_eq()
           << ", ..., " << pn.get_d_eq() << std::endl;
  return output;
}

////////////////////////////////////////////////////////////////////////////////
// ScatteringHabit
////////////////////////////////////////////////////////////////////////////////


ScatteringHabit::ScatteringHabit() {}
ScatteringHabit::~ScatteringHabit() {}

ScatteringHabit::ScatteringHabit(
    const String &name,
    const ArrayOfSingleScatteringData &arts_scat_data,
    const ArrayOfScatteringMetaData &meta_data,
    const Agenda &pnd_agenda,
    const ArrayOfString &pnd_agenda_input)
    : name_(name),
      pnd_agenda_(pnd_agenda),
      pnd_agenda_input_(pnd_agenda_input)
{
  Index n = arts_scat_data.size();
  Index m = meta_data.size();
  std::vector<scattering::Particle> particles{};
  EigenVector d_eq(n), d_max(n), mass(n);

  for (size_t i = 0; i < n; ++i) {
    auto scattering_data = detail::arts_to_scattering(arts_scat_data[i]);
    auto meta = meta_data[i];
    particles.emplace_back(meta.mass,
                           meta.diameter_volume_equ,
                           meta.diameter_max,
                           scattering_data);
  }
  particle_model_ = std::make_shared<scattering::ParticleHabit>(particles);
}

ScatteringHabit::ScatteringHabit(
    const String name,
    const Agenda &pnd_agenda,
    ArrayOfString pnd_agenda_input,
    std::shared_ptr<scattering::ParticleHabit> particle_model)
    : name_(name),
      pnd_agenda_(pnd_agenda),
      pnd_agenda_input_(pnd_agenda_input),
      particle_model_(particle_model) {}

ArrayOfString ScatteringHabit::get_dpnd_data_dx_names(
    ArrayOfRetrievalQuantity jacobian_quantities, bool jacobian_do) const {
  ArrayOfString dpnd_data_dx_names = {};
  if (jacobian_do) {
    for (auto &jq : jacobian_quantities) {
      if (jq.MainTag() == SCATSPECIES_MAINTAG) {
        if (jq.Subtag() == name_) {
          dpnd_data_dx_names.push_back(jq.SubSubtag());
        }
      }
    }
  }
  return dpnd_data_dx_names;
}

Matrix ScatteringHabit::get_agenda_input(Matrix pbp_field,
                                         ArrayOfString pbf_names) const {
    assert(pbp_field.nrows() == pbf_names.size());

  Index n_inputs = pnd_agenda_input_.size();
  ArrayOfIndex row_indices(n_inputs);

  for (Index i = 0; i < n_inputs; ++i) {
    row_indices[i] = find_first(pbf_names, pnd_agenda_input_[i]);
    if (row_indices[i] < 0) {
      ostringstream os;
      os << "Particle habit requires input " << pnd_agenda_input_[i]
         << "\",\nbut this quantity "
         << "could not found in *particle_bulkprop_names*.\n"
         << "(Note that temperature must be written as \"Temperature\")";
      throw runtime_error(os.str());
    }
  }

  Matrix agenda_input(pbp_field.ncols(), n_inputs);
  for (size_t i = 0; i < n_inputs; ++i) {
      agenda_input(joker, i) = pbp_field(row_indices[i], joker);
  }

  return agenda_input;
}

std::shared_ptr<ScatteringSpeciesImpl> ScatteringHabit::prepare_scattering_data(ScatteringPropertiesSpec specs) const {

    scattering::ParticleHabit formatted = particle_model_->set_stokes_dim(specs.n_stokes);

    if (specs.frame == ReferenceFrame::Lab) {
        Index n = specs.lat_scat.nelem();
        EigenVector lon_scat(n);
        Numeric dx = 2.0 * M_PI / static_cast<Numeric>(n - 1);
        lon_scat[0] = 0.0;
        for (Index i = 1; i < n; ++i) {
            lon_scat[i] = lon_scat[i - 1] + dx;
        }
        formatted = formatted.to_lab_frame(to_eigen(specs.lat_inc),
                                           lon_scat,
                                           to_eigen(specs.lat_scat),
                                           specs.n_stokes);
    }

    if (specs.format == Format::Spectral) {
        // Data must be regridded to ensure conformity with grids expected
        // by SHT transform.
        formatted = formatted.regrid().to_spectral(specs.l_max, specs.m_max);
    } else {
        formatted = formatted.downsample_scattering_angles(to_eigen(specs.lon_scat),
                                                           to_eigen(specs.lat_scat));
        formatted = formatted.to_gridded(to_eigen(specs.lon_inc),
                                         to_eigen(specs.lat_inc),
                                         to_eigen(specs.lon_scat),
                                         to_eigen(specs.lat_scat));

    }

    auto new_model = std::make_shared<scattering::ParticleHabit>(formatted.interpolate_frequency(to_eigen(specs.f_grid)));
    auto result = std::make_shared<ScatteringHabit>(name_, pnd_agenda_, pnd_agenda_input_, new_model);
    result->set_phase_function_norm(specs.phase_function_norm);
    return result;
}

BulkScatteringProperties ScatteringHabit::calculate_bulk_properties(
    Workspace &ws,
    ConstMatrixView pbp_field,
    const ArrayOfString &pbf_names,
    ConstVectorView temperature,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
    bool jacobian_do) const {
  Matrix pnd_data;
  Tensor3 dpnd_data_dx;

  Matrix input = get_agenda_input(pbp_field, pbf_names);
  ArrayOfString dpnd_data_dx_names =
      get_dpnd_data_dx_names(jacobian_quantities, jacobian_do);

  pnd_agendaExecute(ws,
                    pnd_data,
                    dpnd_data_dx,
                    *this,
                    temperature,
                    input,
                    pnd_agenda_input_,
                    dpnd_data_dx_names,
                    pnd_agenda_);

  auto n_levels = pnd_data.nrows();
  Array<scattering::SingleScatteringData> bulk_properties(n_levels);
  for (Index i = 0; i < n_levels; ++i) {
      EigenVector number_densities = to_eigen(pnd_data(i, joker));
    bulk_properties[i] = particle_model_->calculate_bulk_properties(
        temperature[i], number_densities, phase_function_norm_);
  }

  return BulkScatteringProperties(bulk_properties);
}

std::ostream &operator<<(std::ostream &output, const ScatteringHabit &habit) {
  output << "Scattering habit: " << habit.name_ << std::endl;
  output << "\t Particle d_eq: " << habit.get_particle_d_eq() << std::endl;
  output << "\t Particle d_max: " << habit.get_particle_d_max() << std::endl;
  output << std::endl;
  return output;
}
