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

  \brief Workspace methods for the manipulation of scatter_species.

*/
#include "arts.h"
#include "scattering_habit.h"
#include "scattering_species.h"

void scattering_speciesAddScatteringHabit(
    Workspace &workspace,
    ArrayOfScatteringSpecies &scattering_species,
    const String &name,
    const ArrayOfSingleScatteringData &scattering_data,
    const ArrayOfScatteringMetaData &meta_data,
    const Agenda &pnd_agenda,
    const ArrayOfString &pnd_agenda_input,
    const Verbosity &verbosity) {
    Agenda checked_agenda = pnd_agenda;
    checked_agenda.set_name("pnd_agenda");
    checked_agenda.check(workspace, verbosity);
  auto scattering_habit =
      std::make_shared<ScatteringHabit>(name, scattering_data, meta_data, checked_agenda, pnd_agenda_input);
  scattering_species.emplace_back(scattering_habit);
}
