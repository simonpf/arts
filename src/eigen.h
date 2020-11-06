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
/** \file eigen.h
 *
 * This header file contains bindings connecting the ARTS vector, matrix
 * and tensor types with those of the Eigen library.
 *
 *
 * \author Simon Pfreundschuh
 * \data 2020-09-15
 *
 */

#ifndef __ARTS_EIGEN_H__
#define __ARTS_EIGEN_H__

#include <Eigen/Core>
#include "Eigen/CXX11/Tensor"
#include "matpackV.h"
#include "matpackVII.h"

////////////////////////////////////////////////////////////////////////////////
// Vectors
////////////////////////////////////////////////////////////////////////////////

using EigenVector = Eigen::Matrix<Numeric, 1, -1, Eigen::RowMajor>;
using EigenVectorMap = Eigen::Map<EigenVector>;
using EigenConstVectorMap = Eigen::Map<const EigenVector>;

inline EigenConstVectorMap to_eigen(const Vector &vector) {
  return EigenConstVectorMap{vector.get_c_array(), vector.nelem()};
}

inline Vector to_arts(const EigenVector &vector) {
  auto n = vector.size();
  Vector result(n);
  std::copy(vector.data(), vector.data() + n, result.get_c_array());
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Tensors
////////////////////////////////////////////////////////////////////////////////

template <int rank>
using EigenTensor = Eigen::Tensor<Numeric, rank, Eigen::RowMajor>;
template <int rank>
using EigenConstTensorMap = Eigen::TensorMap<const EigenTensor<rank>>;

inline EigenConstTensorMap<7> to_eigen(const Tensor7 &tensor) {
  std::array<Eigen::Index, 7> dimensions = {tensor.nlibraries(),
                                            tensor.nvitrines(),
                                            tensor.nshelves(),
                                            tensor.nbooks(),
                                            tensor.npages(),
                                            tensor.nrows(),
                                            tensor.ncols()};
  std::cout << "to eigen:" << std::endl;
  std::cout << dimensions[0] << " // ";
  std::cout << dimensions[1] << " // ";
  std::cout << dimensions[2] << " // ";
  std::cout << dimensions[3] << " // ";
  std::cout << dimensions[4] << " // ";
  std::cout << dimensions[5] << " // ";
  std::cout << dimensions[6] << std::endl;
  return EigenConstTensorMap<7>{tensor.get_c_array(), dimensions};
}

inline EigenConstTensorMap<5> const to_eigen(const Tensor5 &tensor) {
  std::array<Eigen::Index, 5> dimensions = {tensor.nshelves(),
                                            tensor.nbooks(),
                                            tensor.npages(),
                                            tensor.nrows(),
                                            tensor.ncols()};
  return EigenConstTensorMap<5>{tensor.get_c_array(), dimensions};
}

#endif
