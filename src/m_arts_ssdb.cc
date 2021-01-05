/* Copyright (C) 2012
   Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   m_arts_ssdb.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-12-04

  \brief  Workspace function to load particles and particle habits from the
          ARTS single scattering database.
*/
#include "scattering_habit.h"
#include "scattering/arts_ssdb.h"

void scattering_particleReadFromARTSSSDB(ScatteringParticle &out,
                                         const String &filename,
                                         const Verbosity &verbosity) {
    auto particle_file = scattering::arts_ssdb::ParticleFile(filename);
    out = particle_file.to_particle();
}

void particle_habitReadFromARTSSSDB(ArrayOfScatteringParticle &out,
                                    const Index &stokes_dim,
                                    const String &path,
                                    const Verbosity &verbosity) {
    auto habit_folder = scattering::arts_ssdb::HabitFolder(path);
    auto n_particles = habit_folder.get_n_particles();
    out.resize(0);
    out.reserve(n_particles);
    for (auto iterator = habit_folder.begin(); iterator != habit_folder.end(); ++iterator) {
        auto particle = (*iterator).to_particle();
        out.push_back(particle.set_stokes_dim(stokes_dim));
    }
}

void ParticleHabitReadFromARTSSSDB() {}
