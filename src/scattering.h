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
#include "scattering/single_scattering_data.h"

#include "jacobian.h"
#include "matpackI.h"
#include "eigen.h"

#ifndef __ARTS_SCATTERING__
#define __ARTS_SCATTERING__

enum class Format {Gridded, Spectral};
enum class ReferenceFrame {ScatteringPlane, Lab};

/** Specification of required scattering properties.
 *
 * This struct holds requirements describing the scattering
 * data required by different scattering solvers.
 */
struct ScatteringPropertiesSpec {
    ScatteringPropertiesSpec(const Vector &f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             Index l_max_,
                             Index m_max=0,
                             Numeric phase_function_norm=1.0);
    ScatteringPropertiesSpec(const Vector &f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             Vector lon_scat_,
                             Vector lat_scat_,
                             Numeric phase_function_norm=1.0);
    ScatteringPropertiesSpec(const Vector &f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             Vector lat_inc_,
                             Vector lon_scat_,
                             Vector lat_scat_,
                             Numeric phase_function_norm=1.0);
    Format format;
    ReferenceFrame frame;
    Index n_stokes = 0;
    Index l_max = 0;
    Index m_max = 0;
    Numeric phase_function_norm = 1.0;
    Vector lon_inc{};
    Vector lat_inc{};
    Vector lon_scat{};
    Vector lat_scat{};
    Vector f_grid{};
    Index n_angs_frame_conversion = 32;
};


/** BulkScatteringProperties
 *
 * The class represents the bulk scattering properties of the atmosphere. It is the
 * result of calculating the scattering properties all point in the atmosphere. The
 * scattering data for each point in the atmospheric grid is represent by a
 * SingleScatteringData object. The data for the whole grid is stored as a 1-dimensional
 * array.
 */
class BulkScatteringProperties {
 public:
  /** Create BulkScatteringProperties
     *
     * @param Array of scattering::SingleScatteringData ojects describing the scatlib properties
     * at each atmospheric position.
     */
  BulkScatteringProperties(Array<scattering::SingleScatteringData> data)
      : data_(data),
        n_freqs_(data[0].get_f_grid().size()),
        stokes_dim_(data[0].get_stokes_dim()) {}

  /** Extracts extinction coefficients from bulk properties.
   *
   * @return Tensor4 with dimensions [n_layers, n_freqs, n_lon_inc, n_lat_inc]
   * containing the extinction coefficient for all layers in the atmosphere,
   * frequencies in f_grid, and incoming azimuth and zenith angles.
   */
  Tensor4 get_extinction_coeff() const;

  /** Extracts extinction matrices from bulk properties.
   *
   * @return Tensor6 with dimensions
   * [n_layers, n_freqs, n_lon_inc, n_lat_inc, n_stokes, n_stokes]
   * containing the extinction matrices for all layers in the atmosphere,
   * frequencies in f_grid, and incoming azimuth and zenith angles.
   */
  Tensor6 get_extinction_matrix() const;
  Tensor6 get_extinction_matrix(Index stokes_dim) const;

  /** Extracts absorption coefficients from bulk properties.
   *
   * @note Note that extracting the extinction matrix
   *
   * @return Tensor5 with dimensions [n_layers, n_freqs, n_lon_inc]
   * containing the absorption coefficient for all layers in the atmosphere,
   * frequencies in f_grid, and incoming azimuth and zenith angles.
   */
  Tensor4 get_absorption_coeff() const;

  /** Extracts absorption vector from bulk properties.
   *
   * @return Tensor5 with dimensions [n_layers, n_freqs, n_lon_inc, stokes]
   * containing the absorption vector for all layers in the atmosphere,
   * frequencies in f_grid, incoming azimuth and zenith angles and the stokes
   * components.
   */
  Tensor5 get_absorption_vector() const;
  Tensor5 get_absorption_vector(Index stokes_dim) const;

  /** Get spectral coefficients.
   *
   * Return spectral components of SHT-transformed scattering matrix computed
   * using orthonormal basis functions.
   *
   * @return Tensor5 containing the spectral SHT coefficients along
   * column, the atmospheric layers along rows and the frquencies along
   * the pages.
   */
  Tensor5 get_spectral_coeffs() const;

  /** Get Legendre coefficients.
   * @return Tensor5 containing Legendre  coefficients along
   * column, the atmospheric layers along rows and the frquencies along
   * the pages.
   */
  Tensor5 get_legendre_coeffs() const;

  /** Get phase function.
   * @return Tensor3 containing the phase functions elements along
   * column, the atmospheric layers along rows and the frquencies along
   * the pages.
   */
  Tensor3 get_phase_matrix() const;

  Tensor7 get_scattering_matrix() const;
  Tensor7 get_scattering_matrix(Index stokes_dim) const;

  BulkScatteringProperties &operator+=(const BulkScatteringProperties &other) {
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
    return *this;
  }

  /** Normalize phase matrix data.
   *
   * Normalizes the scattering-angle integral of phase matrix to the
   * given value.
   *
   * @param norm The desired integral of the phase matrix.
   */
  void normalize(Numeric norm) {
    for (auto &d : data_) {
      d.normalize(norm);
    }
  }

  void downsample_lon_scat(const Vector &lon_scat) {
    auto lon_scat_ptr = std::make_shared<EigenVector>(to_eigen(lon_scat));
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] = data_[i].downsample_lon_scat(lon_scat_ptr);
    }
  }

 private:
  Array<scattering::SingleScatteringData> data_;
  Index n_freqs_;
  Index stokes_dim_;
};

////////////////////////////////////////////////////////////////////////////////
// Abstract interface
////////////////////////////////////////////////////////////////////////////////

/** Abstract interface for scattering species representations.
 *
 * This abstract class defines the generic interface that all 
 *
 */
class ScatteringSpeciesImpl {
public:
    virtual ~ScatteringSpeciesImpl(){};
    ScatteringSpeciesImpl(){};
    ScatteringSpeciesImpl(const ScatteringSpeciesImpl &) = default;

    virtual std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(
        ScatteringPropertiesSpec specs) const = 0;
    virtual BulkScatteringProperties calculate_bulk_properties(
        Workspace &ws,
        ConstMatrixView pbp_field,
        const ArrayOfString &pbf_names,
        ConstVectorView temperature,
        const ArrayOfRetrievalQuantity &jacobian_quantities,
        bool jacobian_do) const = 0;
};

#endif
