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

  \brief Declaration of classes related to the representation of scattering
         habits, i.e. habits of scattering particles and associated PSD>

   This file contains the definitions of the ScatteringParticle,
   ArrayOfScatteringParticles and of the ScatteringHabit class. A scattering
   particle represents a single scattering particle with given size and
   scattering data. A ScatteringHabit is a collection of sucher particles
   together with a particles size distribution (PSD), which defines the
   distribution of the particles in the habit for different atmopheric states.
*/
#include <iostream>

#include "agenda_class.h"
#include "eigen.h"
#include "optproperties.h"
#include "scattering.h"
#include "scattering/particle_habit.h"

#ifndef __ARTS_SCATTERING_HABIT_H__
#define __ARTS_SCATTERING_HABIT_H__

////////////////////////////////////////////////////////////////////////////////
// ScatteringParticle and ArrayOfScatteringParticle
////////////////////////////////////////////////////////////////////////////////

typedef scattering::Particle ScatteringParticle;
typedef Array<ScatteringParticle> ArrayOfScatteringParticle;

std::ostream &operator<<(std::ostream &out, const ScatteringParticle &);
std::ostream &operator<<(std::ostream &out, const ArrayOfScatteringParticle &);

////////////////////////////////////////////////////////////////////////////////
// ScatteringHabit
////////////////////////////////////////////////////////////////////////////////
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
  /** Create ScatteringHabit.
   * @param name The name of the scattering species.
   * @param scat_data The ensemble scattering data describing the particles
   * which make up the habit.
   * @param meta_data The particle meta data corresponding to the particles in
   * scat_data
   * @pnd_agenda The agenda to use to calculate the number densities of each
   * particle.
   * @pnd_agenda_input The names of the properties from pbp_field to use as input
   * for the particle size distributions.
   */
  ScatteringHabit(const String &name,
                  const ArrayOfSingleScatteringData &scat_data,
                  const ArrayOfScatteringMetaData &meta_data,
                  const Agenda &pnd_agenda,
                  const ArrayOfString &pnd_agenda_input);
  /** Create ScatteringHabit.
   * @param name The name of the scattering species.
   * @pnd_agenda The agenda to use to calculate the number densities of each
   * particle.
   * @pnd_agenda_input The names of the properties from pbp_field to use as input
   * for the particle size distributions.
   * @particle_model Particle model object containing the scattering data
   * describing the particles in the habit.
   */
  ScatteringHabit(const String name,
                  const Agenda &pnd_agenda,
                  ArrayOfString pnd_agenda_input,
                  std::shared_ptr<scattering::ParticleHabit> particle_model);

  ScatteringHabit(const ScatteringHabit &) = default;
  ScatteringHabit &operator=(const ScatteringHabit &) = default;
  ScatteringHabit &operator=(ScatteringHabit &&) = default;

  ~ScatteringHabit();

  /// The masses of the particles in the habit.
  Vector get_particle_mass() const {
    return to_arts(particle_model_->get_mass());
  }

  /// The volume equivalent diameters of the particles in the habit.
  Vector get_particle_d_eq() const {
    return to_arts(particle_model_->get_d_eq());
  }

  /// The maximum diameters of the particles in the habit.
  Vector get_particle_d_max() const {
    return to_arts(particle_model_->get_d_max());
  }

  Vector get_particle_area() const {
    return to_arts(particle_model_->get_d_max());
  }

  void set_phase_function_norm(Numeric norm) { phase_function_norm_ = norm; }

  /** Calculate bulk scattering properties for 1D atmosphere.
   *
   * @param pbp_field Matrix view containing the particle bulk properties for
   * all layers in the atmosphere.
   * @param pbp_names The names corresponding to the rows in pbp_field.
   * @param temperature Vector containing the temperatures in the atmosphere.
   * @param jacobian_quantities The quantities for which to compute the jacobian.
   * @param jacobian_do Whether or not to calculate the Jacobian.
   */
  BulkScatteringProperties calculate_bulk_properties(Workspace &ws,
                                                     ConstMatrixView pbp_field,
                                                     const ArrayOfString &pbp_names,
                                                     ConstVectorView temperature,
                                                     const ArrayOfRetrievalQuantity &jacobian_quantities,
                                                     bool jacobian_do) const;


  std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(ScatteringPropertiesSpec specs) const;

  friend std::ostream &operator<<(std::ostream &out, const ScatteringHabit &);

 private:
  String name_;
  Agenda pnd_agenda_;
  ArrayOfString pnd_agenda_input_;
  Numeric phase_function_norm_ = 1.0;
  std::shared_ptr<scattering::ParticleHabit> particle_model_;
};

#endif
