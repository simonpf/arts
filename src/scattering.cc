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
#include "scattering.h"

////////////////////////////////////////////////////////////////////////////////
// ScatteringPropertiesSpec
////////////////////////////////////////////////////////////////////////////////

ScatteringPropertiesSpec::ScatteringPropertiesSpec(const Vector &f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   Index l_max_,
                                                   Index m_max_,
                                                   Numeric phase_function_norm_)
    : format(Format::Spectral),
      n_stokes(stokes_dim),
      l_max(l_max_),
      m_max(m_max_),
      phase_function_norm(phase_function_norm_),
      f_grid(f_grid_) {}

ScatteringPropertiesSpec::ScatteringPropertiesSpec(const Vector &f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   Vector lon_scat_,
                                                   Vector lat_scat_,
                                                   Numeric phase_function_norm_)
    : format(Format::Gridded),
      n_stokes(stokes_dim),
      phase_function_norm(phase_function_norm_),
      lon_inc(1),
      lat_inc(1),
      lon_scat(lon_scat_),
      lat_scat(lat_scat_),
      f_grid(f_grid_) {
  lon_inc[0] = M_PI;
  lat_inc[0] = 0.5 * M_PI;
}

ScatteringPropertiesSpec::ScatteringPropertiesSpec(const Vector &f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   Vector lat_inc_,
                                                   Vector lon_scat_,
                                                   Vector lat_scat_,
                                                   Numeric phase_function_norm_)
    : format(Format::Gridded),
      frame(frame_),
      n_stokes(stokes_dim),
      phase_function_norm(phase_function_norm_),
      lon_inc(1),
      lat_inc(lat_inc_),
      lon_scat(lon_scat_),
      lat_scat(lat_scat_),
      f_grid(f_grid_) {
  lon_inc[0] = M_PI;
}

////////////////////////////////////////////////////////////////////////////////
// BulkScatteringProperties
////////////////////////////////////////////////////////////////////////////////

Matrix BulkScatteringProperties::get_extinction_coefficients() const {
    auto n = data_.size();
    Matrix result(n_freqs_, n);
    for (Index i = 0; i < n; ++i) {
        EigenTensor<6> extinction = data_[i].get_extinction_coeff();
        result(joker, i) =
            VectorView(extinction.data(), Range(0, extinction.size()));
    }
    return result;
}

Tensor5 BulkScatteringProperties::get_extinction_matrix() const {
  auto n = data_.size();
  auto n_angs = data_[0].get_n_lat_inc();
  Tensor5 result(n_freqs_, n, n_angs, stokes_dim_, stokes_dim_, 0.0);
  for (Index i = 0; i < n; ++i) {
    EigenTensor<8> extinction = data_[i].get_extinction_matrix(stokes_dim_);
    for (Index j = 0; j < n_freqs_; ++j) {
      for (Index k = 0; k < n_angs; ++k) {
        for (Index l = 0; l < stokes_dim_; ++l) {
          for (Index m = 0; m < stokes_dim_; ++m) {
            result(j, i, k, l, m) = extinction.coeffRef(j, 0, 0, k, 0, 0, l, m);
          }
        }
      }
    }
  }
  return result;
}

Matrix BulkScatteringProperties::get_absorption_coefficients() const {
    auto n = data_.size();
    Matrix result(n_freqs_, n);
    for (Index i = 0; i < n; ++i) {
        auto absorption = data_[i].get_absorption_coeff();
        result(joker, i) =
            VectorView(absorption.data(), Range(0, absorption.size()));
    }
    return result;
}

Tensor4 BulkScatteringProperties::get_absorption_vector() const {
  auto n = data_.size();
  auto n_angs = data_[0].get_n_lat_inc();
  Tensor4 result(n_freqs_, n, n_angs, stokes_dim_);
  for (Index i = 0; i < n; ++i) {
    EigenTensor<7> absorption = data_[i].get_absorption_vector(stokes_dim_);
    for (Index j = 0; j < n_freqs_; ++j) {
      for (Index k = 0; k < n_angs; ++k) {
        for (Index l = 0; l < stokes_dim_; ++l) {
          result(j, i, k, l) = absorption.coeffRef(j, 0, 0, k, 0, 0, l);
        }
      }
    }
  }
  return result;
}

  Tensor3 BulkScatteringProperties::get_spectral_coefficients() const {
    auto n = data_.size();
    Index n_coeffs = data_[0].get_phase_function_spectral().dimension(4);
    Tensor3 result(n_freqs_, n, n_coeffs);
    for (Index i = 0; i < n; ++i) {
      EigenTensor<5> phase_matrix =
          data_[i].get_phase_function_spectral().real();
      result(joker, i, joker) = MatrixView(phase_matrix.data(),
                                           Range(0, n_freqs_, n_coeffs),
                                           Range(0, n_coeffs, 1));
    }
    return result;
  }

  Tensor3 BulkScatteringProperties::get_legendre_coefficients() const {
    auto n = data_.size();
    Index n_coeffs = data_[0].get_phase_function_spectral().dimension(4);
    Tensor3 result(n_freqs_, n, n_coeffs);
    for (Index i = 0; i < n; ++i) {
        EigenTensor<5> phase_matrix = data_[i].get_phase_function_spectral().real();
        result(joker, i, joker) = MatrixView(phase_matrix.data(),
                                             Range(0, n_freqs_, n_coeffs),
                                             Range(0, n_coeffs, 1));
    }
    // Need to go from SHT coefficients to Legendre coefficients.
    for (Index i = 0; i < n_coeffs; ++i) {
        result(joker, joker, i) *= sqrt(4.0 * M_PI / (2.0 * i + 1.0));
    }
    return result;
  }

  Tensor6 BulkScatteringProperties::get_scattering_matrix() const {
    auto n_layers = data_.size();
    auto n_lat_inc = data_[0].get_n_lat_inc();
    auto n_lat_scat = data_[0].get_n_lat_scat();
    Tensor6 result(
        n_freqs_, n_layers, n_lat_inc, n_lat_scat, stokes_dim_, stokes_dim_, 0.0);
    for (Index i = 0; i < n_layers; ++i) {
      EigenTensor<8> extinction = data_[i].get_scattering_matrix(stokes_dim_);
      for (Index j = 0; j < n_freqs_; ++j) {
        for (Index k = 0; k < n_lat_inc; ++k) {
          for (Index l = 0; l < n_lat_scat; ++l) {
            for (Index m = 0; m < stokes_dim_; ++m) {
              for (Index n = 0; n < stokes_dim_; ++n) {
                  result(j, i, k, l, m, n) =
                      extinction.coeffRef(j, 0, 0, k, 0, l, m, n);
              }
            }
          }
        }
      }
    }
    return result;
  }

  Tensor3 BulkScatteringProperties::get_phase_matrix() const {
    auto n = data_.size();
    Index n_spectral_coeffs_ = data_[0].get_phase_function().dimension(5);
    Tensor3 result(n_freqs_, n, n_spectral_coeffs_);
    for (Index i = 0; i < n; ++i) {
      EigenTensor<6> phase_matrix = data_[i].get_phase_function();
      result(joker, i, joker) =
          MatrixView(phase_matrix.data(),
                     Range(0, n_freqs_, n_spectral_coeffs_),
                     Range(0, n_spectral_coeffs_, 1));
    }
    return result;
  }
