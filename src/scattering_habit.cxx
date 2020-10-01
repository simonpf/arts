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

EigenTensor<7> artsscat_to_scatlib(const Tensor5 &tensor) {
  EigenTensor<5> tensor_eigen = to_eigen(tensor);
  auto result = scatlib::eigen::unsqueeze<5, 6>(tensor_eigen);
  return scatlib::eigen::cycle_dimensions(result);
}

scatlib::SingleScatteringData artsscat_to_scatlib(
    const SingleScatteringData &arts_data) {
  EigenVector f_grid = to_eigen(arts_data.f_grid);
  EigenVector t_grid = to_eigen(arts_data.T_grid);
  EigenVector angle_grid_inc = EigenVector(1);
  angle_grid_inc[0] = 0.0;
  EigenVector lon_inc = to_eigen(arts_data.aa_grid);
  EigenVector lat_inc = to_eigen(arts_data.za_grid);
  auto phase_matrix = to_eigen(arts_data.pha_mat_data);
  auto extinction_matrix = artsscat_to_scatlib(arts_data.ext_mat_data);
  auto absorption_vector = artsscat_to_scatlib(arts_data.abs_vec_data);
  auto backward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);
  auto forward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);

  return scatlib::SingleScatteringData(f_grid,
                                       t_grid,
                                       angle_grid_inc,
                                       angle_grid_inc,
                                       lon_inc,
                                       lat_inc,
                                       phase_matrix,
                                       extinction_matrix,
                                       absorption_vector,
                                       backward_scattering_coeff,
                                       forward_scattering_coeff);
}

}  // namespace detail

ScatteringHabit::ScatteringHabit() {}
ScatteringHabit::ScatteringHabit(
    const ArrayOfSingleScatteringData &arts_scat_data,
    const ArrayOfScatteringMetaData &meta_data,
    const Agenda &psd_agenda)
    : psd_agenda_(&psd_agenda) {
  Index n = arts_scat_data.size();
  std::vector<scatlib::SingleScatteringData> model_data;
  EigenVector d_eq(n), d_max(n), mass(n);

  for (size_t i = 0; i < arts_scat_data.size(); ++i) {
    auto scattering_data = detail::artsscat_to_scatlib(arts_scat_data[i]);
    model_data.push_back(scattering_data);

    auto meta = meta_data[i];
    d_eq[i] = meta.diameter_volume_equ;
    d_max[i] = meta.diameter_max;
    mass[i] = meta.mass;
  }
  particle_model_ = std::make_shared<scatlib::ParticleModel>(d_eq, d_max, mass, model_data);
}

std::ostream & operator<<(std::ostream &output, const ScatteringHabit &habit) {
    std::cout << "Scattering Habit: NICE" << std::endl;
}
