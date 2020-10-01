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

#ifndef __ARTS_SCATTERING__
#define __ARTS_SCATTERING__

enum class Format {Gridded, Spectral};

struct ScatteringPropertiesSpec {
    ScatteringPropertiesSpec(int l_max_, int m_max=0);
    ScatteringPropertiesSpec(Vector lon_scat_, Vector lat_scat_);

    Format format;
    size_t l_max = 0;
    size_t m_max = 0;
    Vector lon_inc{};
    Vector lat_inc{};
    Vector lon_scat{};
    Vector lat_scat{};
};


class ScatteringProperties {
 public:
  ScatteringProperties(Array<scatlib::SingleScatteringData> data)
      : data_(data), n_freqs(data[0].get_f_grid().size()) {}

  Matrix get_extinction_coefficients() {
    auto n = data_.size();
    Matrix result(n_freqs, n);
    for (Index i = 0; i < n; ++i) {
      auto extinction = data_[i].get_extinction_matrix();
      result(joker, i) =
          VectorView(extinction.data(), Range(0, extinction.size()));
    }
    return result;
  }

  Matrix get_absorption_coefficients() {
    auto n = data_.size();
    Matrix result(n_freqs, n);
    for (Index i = 0; i < n; ++i) {
      auto absorption = data_[i].get_absorption_vector();
      result(joker, i) =
          VectorView(absorption.data(), Range(0, absorption.size()));
    }
    return result;
  }

  Matrix get_spectral_coefficients() {
    auto n = data_.size();
    Index n_spectral_coeffs_ =
        data_[0].get_phase_matrix_spectral().dimension(5);
    Matrix result(n_spectral_coeffs_, n);
    for (Index i = 0; i < n; ++i) {
      EigenTensor<5> phase_matrix =
          data_[i].get_phase_matrix_spectral().real().chip<5>(0);
      result(joker, i) =
          VectorView(phase_matrix.data(), Range(0, phase_matrix.size()));
    }
    return result;
  }

  ScatteringProperties &operator+=(const ScatteringProperties &other) {
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
  }

 private:
  Array<scatlib::SingleScatteringData> data_;
  Index n_freqs;
};

class ScatteringSpeciesImpl {
public:
    virtual ~ScatteringSpeciesImpl() {};
    ScatteringSpeciesImpl() {};
    ScatteringSpeciesImpl(const ScatteringSpeciesImpl &) = default;
    virtual Tensor5 get_phase_matrix(Workspace &ws) = 0;

    virtual void prepare_scattering_data(ScatteringPropertiesSpec specs) = 0;
    virtual ScatteringProperties calculate_bulk_properties(Workspace &ws,
                                                           const MatrixView pbp_field,
                                                           const ArrayOfString pbf_names,
                                                           const Vector temperature,
                                                           const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                           bool jacobian_do) const {};

};


#endif
