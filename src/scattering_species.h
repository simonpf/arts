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
#include "array.h"
#include "agenda_class.h"
#include "eigen.h"
#include "optproperties.h"
#include "scatlib/particle_model.h"
#include "scattering.h"
#include "scattering_habit.h"

#ifndef __ARTS_SCATTERING_SPECIES_H__
#define __ARTS_SCATTERING_SPECIES_H__

class ScatteringSpecies {
 public:
  Tensor5 get_phase_matrix(Workspace &ws) {
    return impl_->get_phase_matrix(ws);
  }

  std::shared_ptr<ScatteringSpeciesImpl> impl_;

  friend std::ostream & operator<<(std::ostream &out, const ScatteringSpecies &);

 private:
};

class ArrayOfScatteringSpecies : public Array<ScatteringSpecies> {
public:



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
