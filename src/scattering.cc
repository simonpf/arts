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
  \file   scattering.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-25

  \brief  Implementation of scattering.h
*/

ScatteringPropertiesSpec::ScatteringPropertiesSpec(Vector lon_scat_,
                                                   Vector lat_scat_)
    : format(Format::Gridded),
      lon_inc(Vector(1)),
      lat_inc(Vector(1)),
      lon_scat(lon_scat),
      lat_scat(lat_scat) {
  lon_inc = 0.0;
  lat_inc = 0.0;
}

ScatteringPropertiesSpec::ScatteringPropertiesSpec(Index l_max_,
                                                   Index m_max_)
    : format(Format::Spectral), lon_scat(lon_scat), lat_scat(lat_scat) {}
