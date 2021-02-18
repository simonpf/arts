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
  \file   scattering_species.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-18

  \brief This file defines the ScatteringSpecies and ArrayOfScatteringSpecies
         classes.

 The ScatteringSpecies class provides a generic interface for all objects
 representing scattering data in ARTS. Together with the
 ArrayOfScatteringSpecies class it implements the high-level interface through
 which all other ARTS components should interact with scattering data.
*/
#include "agenda_class.h"
#include "array.h"
#include "eigen.h"
#include "optproperties.h"
#include "scattering/particle_habit.h"
#include "scattering.h"
#include "scattering_habit.h"

#ifndef __ARTS_SCATTERING_SPECIES_H__
#define __ARTS_SCATTERING_SPECIES_H__


////////////////////////////////////////////////////////////////////////////////
// ScatteringSpecies
////////////////////////////////////////////////////////////////////////////////
/** A generic scattering species.
 *
 * The ScatteringSpecies class provides a generic high-level interface for
 * the representation of scattering data in ARTS.
 *
 * A ScatteringSpecies object is a light-weight container that can hold any
 * type of scattering data. The ScatteringSpecies object manages the
 * lifetime of the scattering data. However, copy semantics are shallow,
 * which means that a copied ScatteringSpecies object shares the underlying
 * data with its source.
 */
class ScatteringSpecies {
 public:

  /// Create empty scattering species.
  ScatteringSpecies() {}

  /// Create container for existing scattering data pointer.
  ScatteringSpecies(std::shared_ptr<ScatteringSpeciesImpl> scattering_habit)
      : impl_(scattering_habit) {}

  /** Prepare scattering data to calculate expected properties.
   *
   * This function prepares scattering data for the calculation of specific
   * properties. It returns a new scattering species object which was
   * converted to allow faster calculation of the scattering properties.
   *
   * @param specs ScatteringPropertiesSpec object describing the expected
   * scattering properties.
   * @return A new scattering species object which was prepared to compute
   * the expected scattering properties.
   */
  ScatteringSpecies prepare_scattering_data(ScatteringPropertiesSpec specs) const {
    return ScatteringSpecies(impl_->prepare_scattering_data(specs));
  }

  /** Calculate bulk scattering properties.
   *
   * @param ws The ARTS workspace to used for calculation of the bulk
   * properties.
   * @param pbp_field The particle bulk properties field flattened into a matrix
   * with rows corresponding to the different atmospheric layers.
   * @param pbp_names The names of the particle bulk properties corresponding
   * to the columns in pbp_field
   * @param temperature The atmospheric temperatures flattened into a vector.
   * @jacobian_quantities The quantities for which to compute the jacobian.
   * @jacobian_do The JacobianDo flag.
   * @return A BulkScatteringProperties object containing the bulk scattering
   * properties.
   */
  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    return impl_->calculate_bulk_properties(ws,
                                            pbp_field,
                                            pbp_names,
                                            temperature,
                                            jacobian_quantities,
                                            jacobian_do);
  }

  friend std::ostream& operator<<(std::ostream& out, const ScatteringSpecies&);

 private:
  std::shared_ptr<ScatteringSpeciesImpl> impl_ = nullptr;
};

////////////////////////////////////////////////////////////////////////////////
// ArrayOfScatteringSpecies
////////////////////////////////////////////////////////////////////////////////

/** Array of ScatteringSpecies
 *
 * This is simply a container class for multiple scattering species. It
 * essentially extends the functions required to calculate bulk scattering
 * properties to an array containing multiple ScatteringSpecies objects.
 */
class ArrayOfScatteringSpecies : public Array<ScatteringSpecies> {
 public:

  ArrayOfScatteringSpecies() = default;
  ArrayOfScatteringSpecies(size_t n) : Array<ScatteringSpecies>(n) {}

  /** Prepare scattering data to calculate expected properties.
   *
   * This function prepares scattering data for the calculation of specific
   * properties. It returns a new scattering species object which was
   * converted to allow faster calculation of the scattering properties.
   *
   * @param specs ScatteringPropertiesSpec object describing the expected
   * scattering properties.
   * @return A new scattering species object which was prepared to compute
   * the expected scattering properties.
   */
  ArrayOfScatteringSpecies prepare_scattering_data(
      ScatteringPropertiesSpec specs) const {
    ArrayOfScatteringSpecies result(size());
    for (size_t i = 0; i < size(); ++i) {
      result[i] = this->operator[](i).prepare_scattering_data(specs);
    }
    result.prepared_ = true;
    result.phase_function_norm_ = specs.phase_function_norm;
    return result;
  }

  /** Calculate bulk scattering properties.
   *
   *
   * @param ws The ARTS workspace to used for calculation of the bulk
   * properties.
   * @param pbp_field The particle bulk properties field flattened into a matrix
   * with rows corresponding to the different atmospheric layers.
   * @param pbp_names The names of the particle bulk properties corresponding
   * to the columns in pbp_field
   * @param temperature The atmospheric temperatures flattened into a vector.
   * @jacobian_quantities The quantities for which to compute the jacobian.
   * @jacobian_do The JacobianDo flag.
   * @return A BulkScatteringProperties object containing the bulk scattering
   * properties.
   */
  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbf_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    if (!prepared_) {
      std::runtime_error(
          "The scattering species must be prepared using "
          "'prepare_scattering' data before the bulk "
          " properties can be computed.");
    }
    auto result =
        this->operator[](0).calculate_bulk_properties(ws,
                                                      pbp_field,
                                                      pbf_names,
                                                      temperature,
                                                      jacobian_quantities,
                                                      jacobian_do);
    for (size_t i = 1; i < this->size(); ++i) {
      result +=
          this->operator[](i).calculate_bulk_properties(ws,
                                                        pbp_field,
                                                        pbf_names,
                                                        temperature,
                                                        jacobian_quantities,
                                                        jacobian_do);
    }
    std::cout << "phase function norm: " << phase_function_norm_ << std::endl;
    result.normalize(phase_function_norm_);
    return result;
  }

private:

  bool prepared_ = false;
  Numeric phase_function_norm_ = 4.0 * M_PI;

};

void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ArrayOfScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity);

void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity);

void Select(  // WS Generic Output:
    ArrayOfScatteringSpecies& needles,
    // WS Generic Input:
    const ArrayOfScatteringSpecies& haystack,
    const ArrayOfIndex& needleind,
    const Verbosity& verbosity);

#endif
