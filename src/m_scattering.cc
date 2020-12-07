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
  \file   scattering.cxx
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-11-27

  \brief Workspace methods to calculate scattering properties.
*/
#include "arts.h"
#include "scattering_species.h"

void scattering_speciesCalcBulkAbsorptionCoeff(Workspace &ws,
                                               Tensor6 &out,
                                               const Vector &f_grid,
                                               const Tensor4 &pbp_field,
                                               const ArrayOfString &pbp_names,
                                               const ArrayOfScatteringSpecies &scattering_species,
                                               const Index &stokes_dim,
                                               const Tensor3 &t_field,
                                               const Vector &aa_grid,
                                               const Vector &za_grid,
                                               const Verbosity &verbosity) {

    Index n_freq = f_grid.nelem();
    Index n_lon_inc = aa_grid.nelem();
    Index n_lat_inc = za_grid.nelem();
    Index n_props = pbp_names.size();
    Index n_p = t_field.npages();
    Index n_lat = t_field.nrows();
    Index n_lon = t_field.ncols();
    Index n_layers = n_lon * n_lat * n_p;

    //
    // Prepare scattering data.
    //

    ScatteringPropertiesSpec scattering_specs(f_grid,
                                              ReferenceFrame::Lab,
                                              stokes_dim,
                                              za_grid,
                                              aa_grid,
                                              za_grid,
                                              1.0);
    auto scattering_species_prepd = scattering_species.prepare_scattering_data(scattering_specs);

    //
    // Flatten atmospheric fields.
    //

    Matrix pbp_field_flat{static_cast<Index>(pbp_names.size()), n_layers};
    Vector temperature = Vector(n_layers);
    Index index = 0;
    for (Index i_p = 0; i_p < n_p; ++i_p) {
        for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
            for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
                pbp_field_flat(joker, index) = pbp_field(joker, i_p, i_lat, i_lon);
                temperature[index] = t_field(i_p, i_lat, i_lon);
                index += 1;
            }
        }
    }

    //
    // Calculate bulk properties.
    //

    auto bulk_properties = scattering_species_prepd.calculate_bulk_properties(ws,
                                                                             pbp_field_flat,
                                                                             pbp_names,
                                                                             temperature,
                                                                             {},
                                                                             false);
    auto absorption = bulk_properties.get_absorption_coeff();

    //
    // Extract absorption.
    //

    index = 0;
    out = Tensor6(n_p, n_lat, n_lon, n_freq, n_lon_inc, n_lat_inc, 0.0);
    for (Index i_p = 0; i_p < n_p; ++i_p) {
        for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
            for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
                for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
                    for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
                        for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
                            out(i_p,  i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc) =
                                absorption(i_freq, index, i_lon_inc, i_lat_inc);
                        }
                    }
                }
                index += 1;
            }
        }
    }

}

