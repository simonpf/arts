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
  \file   scattering_species.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-18

  \brief The ScatteringSpecies container that holds all scattering substances
         in an atmosphere.


   This file contains the definition of the ScatteringSpecies class which
   represents all scattering agents in the atmosphere.
*/
#include "agenda_class.h"
#include "array.h"
#include "eigen.h"
#include "optproperties.h"
#include "scatlib/particle_model.h"
#include "scattering.h"
#include "scattering_habit.h"

#ifndef __ARTS_SCATTERING_SPECIES_H__
#define __ARTS_SCATTERING_SPECIES_H__

class ScatteringSpecies {
 public:
  Tensor5 get_phase_matrix(Workspace& ws) {
    return impl_->get_phase_matrix(ws);
  }

  ScatteringSpecies() {}

  ScatteringSpecies(std::shared_ptr<ScatteringSpeciesImpl> scattering_habit)
      : impl_(scattering_habit) {}

  ScatteringSpecies prepare_scattering_data(ScatteringPropertiesSpec specs) const {
    return ScatteringSpecies(impl_->prepare_scattering_data(specs));
  }

  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbf_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    return impl_->calculate_bulk_properties(ws,
                                            pbp_field,
                                            pbf_names,
                                            temperature,
                                            jacobian_quantities,
                                            jacobian_do);
  }

  friend std::ostream& operator<<(std::ostream& out, const ScatteringSpecies&);

 private:
  std::shared_ptr<ScatteringSpeciesImpl> impl_ = nullptr;
};

class ArrayOfScatteringSpecies : public Array<ScatteringSpecies> {
 public:
  ArrayOfScatteringSpecies() = default;
  ArrayOfScatteringSpecies(size_t n) : Array<ScatteringSpecies>(n) {}

  ArrayOfScatteringSpecies prepare_scattering_data(ScatteringPropertiesSpec specs) const {
    ArrayOfScatteringSpecies result(size());
    for (size_t i = 0; i < size(); ++i) {
        result[i] = this->operator[](i).prepare_scattering_data(specs);
    }
    return result;
  }

  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbf_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    auto result =
        this->operator[](0).calculate_bulk_properties(ws,
                                                      pbp_field,
                                                      pbf_names,
                                                      temperature,
                                                      jacobian_quantities,
                                                      jacobian_do);
    for (size_t i = 1; i < this->size(); ++i) {
      result += this->operator[](i).calculate_bulk_properties(ws,
                                                              pbp_field,
                                                              pbf_names,
                                                              temperature,
                                                              jacobian_quantities,
                                                              jacobian_do);
    }
    return result;
  }
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
