#pragma once

#include "array.h"
#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_muelmat =
    constexpr_smat_data_like<T> and std::remove_cvref_t<T>::size() == 4;

template <typename T>
concept muelmat_convertible = matpack::column_keeper<T> and matpack::row_keeper<T> and matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

struct muelmat final : mat44 {
  constexpr muelmat(Numeric tau = 1.0) noexcept
      : mat44{tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau} {}

  constexpr muelmat(lazy_muelmat auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : mat44(static_cast<mat44>(std::move(a))) {}

  template <muelmat_convertible T>
  muelmat(const T& t) {
    ARTS_USER_ERROR_IF(matpack::column_size(t) != 4, "The matrix must have 4 columns.")
    ARTS_USER_ERROR_IF(matpack::row_size(t) != 4, "The matrix must have 4 rows.")
    for (Index i = 0; i < 4; i++) {
      for (Index j = 0; j < 4; j++) {
        view()(i, j) = matpack::mdvalue(t, {i, j});
      }
    }
  }

  constexpr muelmat(Numeric a, Numeric b, Numeric c, Numeric d, Numeric e,
                    Numeric f, Numeric g, Numeric h, Numeric i, Numeric j,
                    Numeric k, Numeric l, Numeric m, Numeric n, Numeric o,
                    Numeric p) noexcept
      : mat44{a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p} {}
};

//! Addition between muelmat matrices
constexpr auto operator+(const muelmat &a, const muelmat &b) {
  return constexpr_smat_data{a} + constexpr_smat_data{b};
}

//! Addition between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator+(const lazy_muelmat auto &a, const muelmat &b) {
  return a + constexpr_smat_data{b};
}

//! Addition between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator+(const muelmat &a, const lazy_muelmat auto &b) {
  return constexpr_smat_data{a} + b;
}

//! Subtraction between muelmat matrices
constexpr auto operator-(const muelmat &a, const muelmat &b) {
  return constexpr_smat_data{a} - constexpr_smat_data{b};
}

//! Subtraction between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator-(const lazy_muelmat auto &a, const muelmat &b) {
  return a - constexpr_smat_data{b};
}

//! Subtraction between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator-(const muelmat &a, const lazy_muelmat auto &b) {
  return constexpr_smat_data{a} - b;
}

//! Scaling a muelmat matrix
constexpr auto operator*(const Numeric &a, const muelmat &b) {
  return smscl{a, constexpr_smat_data{b}};
}

//! Scaling a muelmat matrix
constexpr auto operator*(const muelmat &a, const Numeric &b) {
  return smscl{b, constexpr_smat_data{a}};
}

//! Scaling a muelmat matrix
constexpr muelmat operator*(const muelmat &a, const muelmat &b) {
  const auto [a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30,
              a31, a32, a33] = a;
  const auto [b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30,
              b31, b32, b33] = b;

  return {a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30,
          a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31,
          a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32,
          a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33,
          a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30,
          a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31,
          a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32,
          a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33,
          a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30,
          a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31,
          a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32,
          a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33,
          a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30,
          a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31,
          a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32,
          a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33};
}

//! Take the average of two muelmat matrices
constexpr auto avg(const muelmat &a, const muelmat &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a muelmat matrix and a non-muelmat matrix
constexpr auto avg(const muelmat &a, const lazy_muelmat auto &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a muelmat matrix and a non-muelmat matrix
constexpr auto avg(const lazy_muelmat auto &a, const muelmat &b) {
  return 0.5 * a + 0.5 * b;
}

//! Get the identity matrix
constexpr muelmat identity_muelmat() { return muelmat{1.0}; }

using muelmat_vector = matpack::matpack_data<muelmat, 1>;
using muelmat_vector_view = matpack::matpack_view<muelmat, 1, false, false>;
using muelmat_vector_const_view = matpack::matpack_view<muelmat, 1, true, false>;

using muelmat_matrix = matpack::matpack_data<muelmat, 2>;
using muelmat_matrix_view = matpack::matpack_view<muelmat, 2, false, false>;
using muelmat_matrix_const_view = matpack::matpack_view<muelmat, 2, true, false>;

using muelmat_tensor3 = matpack::matpack_data<muelmat, 3>;
using muelmat_tensor3_view = matpack::matpack_view<muelmat, 3, false, false>;
using muelmat_tensor3_const_view = matpack::matpack_view<muelmat, 3, true, false>;

muelmat_matrix reverse_cumulative_transmission(const Array<muelmat_vector> &T);

muelmat_matrix forward_cumulative_transmission(const Array<muelmat_vector> &T);
} // namespace rtepack
