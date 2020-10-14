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
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   scattering.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-18

  \brief  Defines the abstract interface for scattering species.
*/
#include "scatlib/single_scattering_data.h"

#include "jacobian.h"
#include "matpackI.h"
#include "eigen.h"

#ifndef __ARTS_SCATTERING__
#define __ARTS_SCATTERING__

enum class Format {Gridded, Spectral};

struct ScatteringPropertiesSpec {
    ScatteringPropertiesSpec(const Vector &f_grid,
                             Index l_max_,
                             Index m_max=0);
    ScatteringPropertiesSpec(const Vector &f_grid,
                             Vector lon_scat_,
                             Vector lat_scat_);

    Format format;
    size_t l_max = 0;
    size_t m_max = 0;
    Vector lon_inc{};
    Vector lat_inc{};
    Vector lon_scat{};
    Vector lat_scat{};

    Vector f_grid{};
};


class ScatteringProperties {
 public:
  ScatteringProperties(Array<scatlib::SingleScatteringData> data)
      : data_(data), n_freqs_(data[0].get_f_grid().size()) {}

  Matrix get_extinction_coefficients() {
    auto n = data_.size();
    Matrix result(n_freqs_, n);
    for (Index i = 0; i < n; ++i) {
      EigenTensor<6> extinction = data_[i].get_extinction_matrix().chip<6>(0);
      result(joker, i) =
          VectorView(extinction.data(), Range(0, extinction.size()));
    }
    return result;
  }

  Matrix get_absorption_coefficients() {
    auto n = data_.size();
    Matrix result(n_freqs_, n);
    for (Index i = 0; i < n; ++i) {
      auto absorption = data_[i].get_absorption_vector();
      result(joker, i) =
          VectorView(absorption.data(), Range(0, absorption.size()));
    }
    return result;
  }

  Tensor3 get_spectral_coefficients() {
    auto n = data_.size();
    Index n_spectral_coeffs_ =
        data_[0].get_phase_matrix_spectral().dimension(5);
    Tensor3 result(n_freqs_, n_spectral_coeffs_, n);
    for (Index i = 0; i < n; ++i) {
      EigenTensor<6> phase_matrix =
          data_[i].get_phase_matrix_spectral().real();
      std::cout << "dim: " << phase_matrix.dimension(0) << std::endl;
      std::cout << "dim: " << phase_matrix.dimension(1) << std::endl;
      std::cout << "dim: " << phase_matrix.dimension(2) << std::endl;
      std::cout << "dim: " << phase_matrix.dimension(3) << std::endl;
      std::cout << "dim: " << phase_matrix.dimension(4) << std::endl;
      std::cout << "dim: " << phase_matrix.dimension(5) << std::endl;
      std::cout << "nspec: " << n_spectral_coeffs_ << std::endl;

    }
    return result;
  }

  ScatteringProperties &operator+=(const ScatteringProperties &other) {
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
    return *this;
  }

 private:
  Array<scatlib::SingleScatteringData> data_;
  Index n_freqs_;
};

class ScatteringSpeciesImpl {
public:
    virtual ~ScatteringSpeciesImpl() {};
    ScatteringSpeciesImpl() {};
    ScatteringSpeciesImpl(const ScatteringSpeciesImpl &) = default;
    virtual Tensor5 get_phase_matrix(Workspace &ws) = 0;

    virtual std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(ScatteringPropertiesSpec specs) const = 0;
    virtual ScatteringProperties calculate_bulk_properties(Workspace &ws,
                                                           const MatrixView pbp_field,
                                                           const ArrayOfString pbf_names,
                                                           const Vector temperature,
                                                           const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                           bool jacobian_do) const = 0;
};


#endif
