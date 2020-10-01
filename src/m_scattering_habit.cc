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
  \date   2020-09-22

  \brief Workspace methods for the manipulation of scattering_habit.
*/
#include "scattering_habit.h"
#include "microphysics.h"

void scattering_habitGetParticleSizes(
    Vector& scat_species_x,
    Numeric& scat_species_a,
    Numeric& scat_species_b,
    const ScatteringHabit& scattering_habit,
    const String& x_unit,
    const Numeric& x_fit_start,
    const Numeric& x_fit_end,
    const Index& do_only_x,
    const Verbosity&) {

  Vector mass = scattering_habit.get_particle_mass();

  if (x_unit == "dveq") {
      scat_species_x = scattering_habit.get_particle_d_eq();
  } else if (x_unit == "dmax") {
      scat_species_x = scattering_habit.get_particle_d_max();
  } else if (x_unit == "area") {
      scat_species_x = scattering_habit.get_particle_area();
  } else if (x_unit == "mass") {
      scat_species_x = mass;
  } else {
      ostringstream os;
      os << "You have selected the x_unit: " << x_unit
         << "while accepted choices are: \"dveq\", \"dmax\", \"mass\" and \"area\"";
      throw runtime_error(os.str());
  }

  if (do_only_x) {
      scat_species_a = -1;
      scat_species_b = -1;
  } else {
      derive_scat_species_a_and_b(scat_species_a,
                                  scat_species_b,
                                  scat_species_x,
                                  mass,
                                  x_fit_start,
                                  x_fit_end);
  }
}
