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
#include <locale>
#include "arts.h"
#include "scattering_species.h"
#include "constants.h"

scattering::LatitudeGridPtr<Numeric> get_quadrature(const String &name,
                                                    const Vector &grid,
                                                    Index degree) {
    if (name == "irregular") {
        return std::make_shared<scattering::IrregularLatitudeGrid<Numeric>>(to_eigen(grid));
    } else if (name == "gauss-legendre") {
        return std::make_shared<scattering::GaussLegendreGrid<Numeric>>(degree);
    } else if (name == "double-gauss") {
        return std::make_shared<scattering::DoubleGaussGrid<Numeric>>(degree);
    } else if (name == "lobatto") {
        return std::make_shared<scattering::LobattoGrid<Numeric>>(degree);
    }
    std::runtime_error("The provided quadrature name is not supported.");
}

ArrayOfScatteringSpecies calculate_bulk_properties_lab(
    const ArrayOfScatteringSpecies &scattering_species,
    const Vector &f_grid,
    const Index &stokes_dim,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    Index quadrature_degree) {

    //
    // Prepare scattering data.
    //

    auto f_grid_ptr = std::make_shared<EigenVector>(to_eigen(f_grid));
    auto aa_grid_ptr = std::make_shared<EigenVector>(to_eigen(aa_grid));
    auto za_grid_ptr = get_quadrature(quadrature_type, za_grid, quadrature_degree);

    ScatteringPropertiesSpec scattering_specs(
        f_grid_ptr, ReferenceFrame::Lab, stokes_dim, za_grid_ptr, aa_grid_ptr,
        za_grid_ptr, 1.0);
    return scattering_species.prepare_scattering_data(scattering_specs);
}

ArrayOfScatteringSpecies calculate_bulk_properties_legendre(
    const ArrayOfScatteringSpecies &scattering_species,
    const Vector &f_grid,
    const Index &stokes_dim,
    const Index &l_max) {

    //
    // Prepare scattering data.
    //

    ScatteringPropertiesSpec scattering_specs(f_grid,
                                              ReferenceFrame::ScatteringPlane,
                                              stokes_dim,
                                              l_max);
    return scattering_species.prepare_scattering_data(scattering_specs);
}


BulkScatteringProperties extract_bulk_properties(
    const ArrayOfScatteringSpecies &scattering_species_prepd,
    Workspace &ws,
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
  Index n_lon_inc = 1;
  Index n_lat_inc = za_grid.nelem();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

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

  return scattering_species_prepd.calculate_bulk_properties(
      ws, pbp_field_flat, pbp_names, temperature, {}, false);
}


void scattering_speciesCalcBulkAbsorptionCoeff(
    Workspace &ws,
    Tensor6 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    const Index &quadrature_degree,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = 1;
  Index n_lat_inc = za_grid.nelem();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

  auto scattering_species_prepd = calculate_bulk_properties_lab(scattering_species, f_grid, stokes_dim, aa_grid,
                                                                za_grid, quadrature_type, quadrature_degree);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto absorption = bulk_properties.get_absorption_coeff();

  //
  // Extract absorption.
  //

  Index index = 0;
  out = Tensor6(n_p, n_lat, n_lon, n_freq, n_lon_inc, n_lat_inc, 0.0);
  for (Index i_p = 0; i_p < n_p; ++i_p) {
    for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
      for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
        for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
          for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
            for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
              out(i_p, i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc) =
                  absorption(i_freq, index, i_lon_inc, i_lat_inc);
            }
          }
        }
        index += 1;
      }
    }
  }
}

void scattering_speciesCalcBulkExtinctionCoeff(
    Workspace &ws,
    Tensor6 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    const Index &quadrature_degree,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = 1;
  Index n_lat_inc = za_grid.nelem();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;


  auto scattering_species_prepd = calculate_bulk_properties_lab(scattering_species, f_grid, stokes_dim, aa_grid,
                                                                za_grid, quadrature_type, quadrature_degree);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto extinction = bulk_properties.get_extinction_coeff();

  //
  // Extract extinction.
  //

  Index index = 0;
  out = Tensor6(n_p, n_lat, n_lon, n_freq, n_lon_inc, n_lat_inc, 0.0);
  for (Index i_p = 0; i_p < n_p; ++i_p) {
    for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
      for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
        for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
          for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
            for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
              out(i_p, i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc) =
                  extinction(i_freq, index, i_lon_inc, i_lat_inc);
            }
          }
        }
        index += 1;
      }
    }
  }
}

void scattering_speciesCalcPhaseFunctionLegendreCoeffs(
    Workspace &ws,
    Tensor7 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const Index &n_coeffs,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = aa_grid.nelem();
  Index n_lat_inc = za_grid.nelem();
  Index n_props = pbp_names.size();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

  auto scattering_species_prepd = calculate_bulk_properties_legendre(scattering_species, f_grid, stokes_dim, n_coeffs);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto coeffs = bulk_properties.get_legendre_coeffs();

  //
  // Extract extinction.
  //

  Index index = 0;
  out = Tensor7(n_p, n_lat, n_lon, n_freq, n_lon_inc, n_lat_inc, n_coeffs, 0.0);
  Index n_coeffs_out = std::min(n_coeffs, coeffs.ncols());
  for (Index i_p = 0; i_p < n_p; ++i_p) {
    for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
      for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
        for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
          for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
            for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
              for (Index i_coeff = 0; i_coeff < n_coeffs_out; ++i_coeff) {
                out(i_p, i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc, i_coeff) =
                    coeffs(i_freq, index, i_lon_inc, i_lat_inc, i_coeff);
              }
            }
          }
        }
        index += 1;
      }
    }
  }
}

void scattering_speciesCalcAbsorptionVector(
    Workspace &ws,
    Tensor7 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    const Index &quadrature_degree,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = aa_grid.nelem();
  Index n_lat_inc = za_grid.nelem();
  Index n_props = pbp_names.size();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

  auto scattering_species_prepd = calculate_bulk_properties_lab(scattering_species, f_grid, stokes_dim,
                                                                aa_grid, za_grid, quadrature_type, quadrature_degree);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto absorption = bulk_properties.get_absorption_vector(stokes_dim);

  //
  // Extract extinction.
  //

  Index index = 0;
  Index n_coeffs = absorption.ncols();
  out = Tensor7(n_p, n_lat, n_lon, n_freq, n_lon_inc, n_lat_inc, n_coeffs, 0.0);
  for (Index i_p = 0; i_p < n_p; ++i_p) {
    for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
      for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
        for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
          for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
            for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
              for (Index i_coeff = 0; i_coeff < n_coeffs; ++i_coeff) {
                out(i_p, i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc, i_coeff) =
                    absorption(i_freq, index, i_lon_inc, i_lat_inc, i_coeff);
              }
            }
          }
        }
        index += 1;
      }
    }
  }
}

void scattering_speciesCalcExtinctionMatrix(
    Workspace &ws,
    Tensor7 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    const Index &quadrature_degree,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = 1;
  Index n_lat_inc = za_grid.nelem();
  Index n_lon_scat = aa_grid.nelem();
  Index n_lat_scat = za_grid.nelem();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

  auto scattering_species_prepd = calculate_bulk_properties_lab(scattering_species, f_grid, stokes_dim,
                                                                aa_grid, za_grid, quadrature_type, quadrature_degree);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto extinction = bulk_properties.get_extinction_matrix(stokes_dim);

  //
  // Extract extinction.
  //

  Index index = 0;
  Index n_coeffs = extinction.ncols();
  out = Tensor7(n_p,
                n_lat,
                n_lon,
                n_freq,
                n_lon_inc,
                n_lat_inc,
                n_coeffs * n_coeffs,
                0.0);
  for (Index i_p = 0; i_p < n_p; ++i_p) {
    for (Index i_lat = 0; i_lat < n_lat; ++i_lat) {
      for (Index i_lon = 0; i_lon < n_lon; ++i_lon) {
        for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
          for (Index i_lon_inc = 0; i_lon_inc < n_lon_inc; ++i_lon_inc) {
            for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
              for (Index i_coeff = 0; i_coeff < n_coeffs * n_coeffs;
                   ++i_coeff) {
                out(i_p, i_lat, i_lon, i_freq, i_lon_inc, i_lat_inc, i_coeff) =
                    extinction(i_freq,
                               index,
                               i_lon_inc,
                               i_lat_inc,
                               i_coeff / n_coeffs,
                               i_coeff % n_coeffs);
              }
            }
          }
        }
        index += 1;
      }
    }
  }
}

void scattering_speciesCalcScatteringMatrix(
    Workspace &ws,
    Tensor7 &out,
    const Vector &f_grid,
    const Tensor4 &pbp_field,
    const ArrayOfString &pbp_names,
    const ArrayOfScatteringSpecies &scattering_species,
    const Index &stokes_dim,
    const Tensor3 &t_field,
    const Vector &aa_grid,
    const Vector &za_grid,
    const String &quadrature_type,
    const Index &quadrature_degree,
    const Verbosity &verbosity) {
  Index n_freq = f_grid.nelem();
  Index n_lon_inc = 1;
  Index n_lat_inc = za_grid.nelem();
  Index n_lon_scat = aa_grid.nelem();
  Index n_lat_scat = za_grid.nelem();
  Index n_p = t_field.npages();
  Index n_lat = t_field.nrows();
  Index n_lon = t_field.ncols();
  Index n_layers = n_lon * n_lat * n_p;

  Vector lon_scat = aa_grid;
  lon_scat *= Conversion::DEG2RAD;
  Vector lat_scat = za_grid;
  lat_scat *= Conversion::DEG2RAD;
  auto scattering_species_prepd = calculate_bulk_properties_lab(scattering_species, f_grid, stokes_dim,
                                                                aa_grid, za_grid, quadrature_type, quadrature_degree);
  auto bulk_properties = extract_bulk_properties(scattering_species_prepd,
                                                 ws,
                                                 f_grid,
                                                 pbp_field,
                                                 pbp_names,
                                                 scattering_species,
                                                 stokes_dim,
                                                 t_field,
                                                 aa_grid,
                                                 za_grid,
                                                 verbosity);
  auto z = bulk_properties.get_scattering_matrix(stokes_dim);

  //
  // Extract extinction.
  //

  Index index = 0;
  Index n_coeffs = z.ncols();
  out = Tensor7(n_layers,
                n_freq,
                n_lat_inc,
                n_lon_scat,
                n_lat_scat,
                n_coeffs,
                n_coeffs,
                0.0);
  for (Index i_l = 0; i_l < n_layers; ++i_l) {
    for (Index i_freq = 0; i_freq < n_freq; ++i_freq) {
      for (Index i_lat_inc = 0; i_lat_inc < n_lat_inc; ++i_lat_inc) {
        for (Index i_lon_scat = 0; i_lon_scat < n_lon_scat; ++i_lon_scat) {
          for (Index i_lat_scat = 0; i_lat_scat < n_lat_scat; ++i_lat_scat) {
            for (Index i_coeff = 0; i_coeff < n_coeffs; ++i_coeff) {
              for (Index j_coeff = 0; j_coeff < n_coeffs; ++j_coeff) {
                out(i_l,
                    i_freq,
                    i_lat_inc,
                    i_lon_scat,
                    i_lat_scat,
                    i_coeff,
                    j_coeff) = z(i_freq,
                                 i_l,
                                 i_lat_inc,
                                 i_lon_scat,
                                 i_lat_scat,
                                 i_coeff,
                                 j_coeff);
              }
            }
          }
        }
      }
    }
  }
}
