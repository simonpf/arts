/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmer.se>

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
  \file   scattering_habit.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-15

  \brief  Representation of scattering particle species with given
          functional size distribution.

   This file contains the definition of the ScatteringHabit class which
   represents a scattering species which is described by an explicit scattering
   model consisting of scatter properties for different particle sizes and
   a particle size distribution.
*/
#include <iostream>

#include "agenda_class.h"
#include "eigen.h"
#include "optproperties.h"
#include "scattering.h"
#include "matpackV.h"
#include "scatlib/particle_model.h"

#ifndef __ARTS_SCATTERING_HABIT_H__
#define __ARTS_SCATTERING_HABIT_H__

namespace detail {

/** Extract backward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scatlib format.
 * @return Eigen tensor containing the backscattering coefficient in
 *     scatlib-compatible format.
 */
EigenTensor<7> extract_backward_scattering_coeff(
    const EigenTensor<7> &phase_matrix);

/** Extract forward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scatlib format.
 * @return Eigen tensor containing the forward scattering coefficient in
 *     scatlib-compatible format.
 */
EigenTensor<7> extract_forward_scattering_coeff(
    const EigenTensor<7> &phase_matrix);

/** Convert ARTS to scatlib format.
 *
 * This function converts ARTS phase matrix data to scatlib format.
 *
 * @param tensor ARTS phase matrix data.
 * @return Eigen tensor containing the phase matrix in scatlib-compatible format.
 */
EigenTensor<7> arts_to_scatlib(const Tensor7 &tensor);

/** Convert ARTS to scatlib format.
 *
 * This function converts a ARTS extinction matrix or absorption
 * vector data to scatlib format
 *
 * @param tensor The input data to convert.
 * @return Eigen tensor containing the data in scatlib-compatible format.
 */
EigenTensor<7> arts_to_scatlib(const Tensor5 &tensor);

scatlib::SingleScatteringData arts_to_scatlib(
    const SingleScatteringData &arts_data);


}  // namespace detail

/** ScatteringHabit
 *
 * A scattering habit describes an atmospheric scattering species through  an
 * explicit particle model of scattering properties at given particle sizes and
 * a parametrized PSD that describes the distribution of those sizes throughout
 * the atmosphere.
 *
 * This class implements the abstract interface for scattering species (defined
 * by ScatteringSpeciesImpl) and instances of this class are used as a elements
 * of the *scattering_species* WSM containing the atmosphere's scattering
 * species.
 *
 */
class ScatteringHabit : public ScatteringSpeciesImpl {

    // Extracts agenda input from particle bulkprop field.
    Matrix get_agenda_input(Matrix pbp_field,
                            ArrayOfString pbf_names) const;

    // Get names of agenda input for which to calculate Jacobians.
    ArrayOfString get_dpnd_data_dx_names(ArrayOfRetrievalQuantity jacobian_quantities,
                                         bool jacobian_do) const;

 public:
  ScatteringHabit();
  /** Create ScatteringHabit
   *
   * @param name The name of the scattering species.
   * @param arts_scat_data The ensemble scattering data describing the particles
   * that make up the habit.
   * @param The particle meta data corresponding to the particles in arts_scat_data
   * @
   *

   */

  ScatteringHabit(const String &name,
                  const ArrayOfSingleScatteringData &arts_scat_data,
                  const ArrayOfScatteringMetaData &meta_data,
                  const Agenda &pnd_agenda,
                  const ArrayOfString &pnd_agenda_input);
  ScatteringHabit(const String name,
                  const Agenda &pnd_agenda,
                  ArrayOfString pnd_agenda_input,
                  std::shared_ptr<scatlib::ParticleModel> particle_model);

  ScatteringHabit(const ScatteringHabit &) = default;
  ScatteringHabit &operator=(const ScatteringHabit &) = default;
  ScatteringHabit &operator=(ScatteringHabit &&) = default;

  ~ScatteringHabit();

  Tensor5 get_phase_matrix(Workspace &ws) { return Tensor5(); }

  Vector get_particle_mass() const {
    return to_arts(particle_model_->get_mass());
  }

  Vector get_particle_d_eq() const {
    return to_arts(particle_model_->get_d_eq());
  }

  Vector get_particle_d_max() const {
    return to_arts(particle_model_->get_d_max());
  }

  Vector get_particle_area() const {
    return to_arts(particle_model_->get_d_max());
  }

  void set_phase_function_norm(Numeric norm) { phase_function_norm_ = norm; }

  BulkScatteringProperties calculate_bulk_properties(Workspace &ws,
                                                     const MatrixView pbp_field,
                                                     const ArrayOfString pbf_names,
                                                     const Vector temperature,
                                                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                     bool jacobian_do) const;

  std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(ScatteringPropertiesSpec specs) const;

  friend std::ostream &operator<<(std::ostream &out, const ScatteringHabit &);

 private:
  String name_;
  Agenda pnd_agenda_;
  ArrayOfString pnd_agenda_input_;
  Numeric phase_function_norm_ = 1.0;
  std::shared_ptr<scatlib::ParticleModel> particle_model_;
};


#endif
