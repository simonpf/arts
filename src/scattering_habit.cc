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
  \file   scattering_habit.cxx
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-15

  \brief  Implementation of scattering_habit.h
*/
#include "scattering_habit.h"
#include "auto_md.h"

const String SCATSPECIES_MAINTAG = "Scattering species";

namespace detail {

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

EigenTensor<7> arts_to_scatlib(const Tensor7 &tensor) {
    EigenTensor<7> tensor_eigen = to_eigen(tensor);
    std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 5, 4, 3, 2, 6};
    return tensor_eigen.shuffle(shuffle_dimensions);
}

EigenTensor<7> arts_to_scatlib(const Tensor5 &tensor) {
  EigenTensor<5> tensor_eigen = to_eigen(tensor);
  auto result = scatlib::eigen::unsqueeze<4, 5>(tensor_eigen);
  std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 3, 2, 4, 5, 6};
  return result.shuffle(shuffle_dimensions);
}

scatlib::SingleScatteringData artsscat_to_scatlib(
    const SingleScatteringData &arts_data) {
  EigenVector f_grid = to_eigen(arts_data.f_grid);
  EigenVector t_grid = to_eigen(arts_data.T_grid);
  EigenVector angle_grid_inc = EigenVector(1);
  angle_grid_inc[0] = 0.0;
  EigenVector lon_scat = to_eigen(arts_data.aa_grid);
  if (lon_scat.size() == 0) {
      lon_scat = angle_grid_inc;
  }
  EigenVector lat_scat = to_eigen(arts_data.za_grid);
  auto phase_matrix = arts_to_scatlib(arts_data.pha_mat_data);
  auto extinction_matrix = arts_to_scatlib(arts_data.ext_mat_data);
  auto absorption_vector = arts_to_scatlib(arts_data.abs_vec_data);
  auto backward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);
  auto forward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);

  auto dimensions = phase_matrix.dimensions();

  assert(phase_matrix.dimension(0) == f_grid.size());
  assert(phase_matrix.dimension(1) == t_grid.size());
  assert(phase_matrix.dimension(4) == lon_scat.size());
  assert(phase_matrix.dimension(5) == lat_scat.size());

  return scatlib::SingleScatteringData(f_grid,
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

//ScatteringHabit::ScatteringHabit() {}
//ScatteringHabit::~ScatteringHabit() {}
ScatteringHabit::ScatteringHabit(
    const String &name,
    const ArrayOfSingleScatteringData &arts_scat_data,
    const ArrayOfScatteringMetaData &meta_data,
    const Agenda &pnd_agenda,
    const ArrayOfString &pnd_agenda_input)
    : name_(name),
      pnd_agenda_(&pnd_agenda),
      pnd_agenda_input_(pnd_agenda_input)
{
  Index n = arts_scat_data.size();
  Index m = meta_data.size();
  std::cout << n << " / " << m << std::endl;
  std::vector<scatlib::SingleScatteringData> model_data;
  EigenVector d_eq(n), d_max(n), mass(n);

  for (size_t i = 0; i < n; ++i) {
    auto scattering_data = detail::artsscat_to_scatlib(arts_scat_data[i]);
    model_data.push_back(scattering_data);

    auto meta = meta_data[i];
    d_eq[i] = meta.diameter_volume_equ;
    d_max[i] = meta.diameter_max;
    mass[i] = meta.mass;
  }
  particle_model_ = std::make_shared<scatlib::ParticleModel>(d_eq, d_max, mass, model_data);
}

ArrayOfString ScatteringHabit::get_dpnd_data_dx_names(
    ArrayOfRetrievalQuantity jacobian_quantities, bool jacobian_do) {
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
                                         ArrayOfString pbf_names) {
  assert(pbp_field.ncols() == pbf_names.size());

  Index n_inputs = pnd_agenda_input_.size();
  ArrayOfIndex column_indices(n_inputs);

  for (Index i = 0; i < n_inputs; ++i) {
    column_indices[i] = find_first(pbf_names, pnd_agenda_input_[i]);
    if (column_indices[i] < 0) {
      ostringstream os;
      os << "Particle habit requires input " << pnd_agenda_input_[i]
         << "\",\nbut this quantity "
         << "could not found in *particle_bulkprop_names*.\n"
         << "(Note that temperature must be written as \"Temperature\")";
      throw runtime_error(os.str());
    }
  }

  Matrix agenda_input(pbp_field.nrows(), n_inputs);
  for (size_t i = 0; i < n_inputs; ++i) {
      agenda_input(joker, i) = pbp_field(joker, column_indices[i]);
  }

  return agenda_input;
}

//ScatteringProperties ScatteringHabit::calculate_bulk_properties(Workspace &ws,
//const MatrixView pbp_field,
//const ArrayOfString pbf_names,
//const Vector temperature,
//const ArrayOfRetrievalQuantity &jacobian_quantities,
//bool jacobian_do) {
//
//Matrix pnd_data;
//Tensor3 dpnd_data_dx;
//
//
//Matrix input = get_agenda_input(pbp_field, pbf_names);
//ArrayOfString dpnd_data_dx_names = get_dpnd_data_dx_names(jacobian_quantities,
//jacobian_do);
//
//pnd_agendaExecute(ws,
//pnd_data,
//dpnd_data_dx,
//temperature,
//input,
//pnd_agenda_input_,
//dpnd_data_dx_names,
                      //*pnd_agenda_);
//
//auto n_levels = pnd_data.nrows();
//Array<scatlib::SingleScatteringData> bulk_properties(n_levels);
//for (Index i = 0; i < n_levels; ++i) {
//bulk_properties[i] = particle_model_->calculate_bulk_properties(temperature[i],
//to_eigen(pnd_data(i, joker)));
//}

    return ScatteringProperties(bulk_properties);

}

std::ostream & operator<<(std::ostream &output, const ScatteringHabit &habit) {
    std::cout << "Scattering Habit: NICE" << std::endl;
}
