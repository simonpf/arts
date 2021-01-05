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
#include <filesystem>

#include "scattering_habit.h"
#include "microphysics.h"
#include "xml_io.h"

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
  std::cout << scattering_habit << std::endl;

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

void scattering_particleReadFromLegacyFormat(
    ScatteringParticle& out,
    const String& single_scattering_data_file,
    const String& scattering_meta_data_file,
    const String& name,
    const Verbosity& verbosity) {
  SingleScatteringData scattering_data;
  xml_read_from_file(single_scattering_data_file, scattering_data, verbosity);
  ScatteringMetaData meta_data;
  xml_read_from_file(scattering_meta_data_file, meta_data, verbosity);

  scattering::ParticleProperties properties{
      name,
      meta_data.source,
      meta_data.refr_index,
      meta_data.mass,
      meta_data.diameter_volume_equ,
      meta_data.diameter_max,
      meta_data.diameter_area_equ_aerodynamical};
  out = ScatteringParticle(properties,
                           detail::from_legacy_format(scattering_data));
}

void scattering_particlesReadFromLegacyFormat(
    ArrayOfScatteringParticle& out,
    const String& array_of_single_scattering_data_file,
    const String& array_of_scattering_meta_data_file,
    const String& name,
    const Verbosity& verbosity) {
  ArrayOfSingleScatteringData scattering_data;
  xml_read_from_file(
      array_of_single_scattering_data_file, scattering_data, verbosity);
  ArrayOfScatteringMetaData meta_data;
  xml_read_from_file(array_of_scattering_meta_data_file, meta_data, verbosity);

  std::string particle_name = name;
  if (particle_name == "") {
    particle_name =
        std::filesystem::path(
            static_cast<std::string>(array_of_single_scattering_data_file))
            .stem();
  }
  out = ArrayOfScatteringParticle();
  for (Index i = 0; i < scattering_data.nelem(); ++i) {
    auto m = meta_data[i];
    scattering::ParticleProperties properties{
        particle_name,
        m.source,
        m.refr_index,
        m.mass,
        m.diameter_volume_equ,
        m.diameter_max,
        m.diameter_area_equ_aerodynamical};
    out.push_back(ScatteringParticle(
        properties, detail::from_legacy_format(scattering_data[i])));
  }
}
