// Orthonormal Dubiner modal basis on the reference triangle.
//
// The reference triangle has vertices
//   v0 = (-1, -1),  v1 = ( 1, -1),  v2 = (-1,  1),
// i.e. the set {(xi, eta) : xi >= -1, eta >= -1, xi + eta <= 0}, with area 2.
//
// The basis uses the Dubiner collapsed-coordinate construction. Map each point
// to the reference square (a, b) in [-1, 1] x [-1, 1] by
//   b = eta,
//   a = (2 + 2 xi) / (1 - eta) - 1,   with a = -1 at the apex eta = 1.
// In these coordinates the orthonormal mode (i, j), i, j >= 0, i + j <= p, is
//   phi_ij(xi, eta) = N_ij * P_i^(0,0)(a) * ((1 - b) / 2)^i * P_j^(2i+1,0)(b),
//   N_ij = sqrt((2 i + 1) (i + j + 1) / 2),
// where P_n^(alpha,beta) is the standard Jacobi polynomial. With this
// normalization
//   integral_T phi_ij phi_kl dA = delta_{ik} delta_{jl}.
// Mode zero is phi_00 = 1 / sqrt(2), the only mode with nonzero mean, so it
// carries the cell average (matching the orthonormal Legendre line basis).
//
// Modes are ordered by total degree n = i + j ascending, then by the first
// reference-coordinate exponent i ascending: within degree n the order is
// (0, n), (1, n-1), ..., (n, 0). The count through degree p is
// (p + 1)(p + 2) / 2.
//
// These evaluators are header-only and allocation-free.

#pragma once

#include <cmath>

#include "wasl/basis/jacobi.hpp"

namespace wasl {
namespace detail {

// Collapsed coordinate a(xi, eta); the apex eta = 1 maps to a = -1.
inline double collapsed_a(double xi, double eta) noexcept {
  return (eta < 1.0) ? (2.0 + 2.0 * xi) / (1.0 - eta) - 1.0 : -1.0;
}

inline double dubiner_norm(int i, int n) noexcept {
  return std::sqrt((2.0 * i + 1.0) * (n + 1.0) / 2.0);
}

} // namespace detail

// Evaluate orthonormal Dubiner modes 0..(triangle_mode_count(p)-1) at (xi, eta)
// into values, which must hold at least triangle_mode_count(p) doubles.
inline void dubiner_basis(double xi, double eta, int p, double *values) noexcept {
  const double b = eta;
  const double a = detail::collapsed_a(xi, eta);
  const double half_one_minus_b = 0.5 * (1.0 - b);

  int index = 0;
  for (int n = 0; n <= p; ++n) {
    for (int i = 0; i <= n; ++i) {
      const int j = n - i;
      const double pi = jacobi_p(a, i, 0.0, 0.0);
      const double warp = std::pow(half_one_minus_b, static_cast<double>(i));
      const double qj = jacobi_p(b, j, 2.0 * i + 1.0, 0.0);
      values[index++] = detail::dubiner_norm(i, n) * pi * warp * qj;
    }
  }
}

// Evaluate orthonormal Dubiner gradients (d/dxi, d/deta) at (xi, eta). Each
// output must hold at least triangle_mode_count(p) doubles. Gradients are
// defined for interior points (eta < 1); the apex is a collapsed-coordinate
// singularity.
inline void dubiner_basis_gradient(double xi, double eta, int p, double *gradients_xi,
                                   double *gradients_eta) noexcept {
  const double b = eta;
  const double a = detail::collapsed_a(xi, eta);
  const double half_one_minus_b = 0.5 * (1.0 - b);
  const double inverse_one_minus_b = 1.0 / (1.0 - b);

  // Collapsed-coordinate chain-rule factors.
  const double da_dxi = 2.0 * inverse_one_minus_b;
  const double da_deta = (a + 1.0) * inverse_one_minus_b;

  int index = 0;
  for (int n = 0; n <= p; ++n) {
    for (int i = 0; i <= n; ++i) {
      const int j = n - i;
      const double norm = detail::dubiner_norm(i, n);
      const double pi = jacobi_p(a, i, 0.0, 0.0);
      const double warp = std::pow(half_one_minus_b, static_cast<double>(i));
      const double qj = jacobi_p(b, j, 2.0 * i + 1.0, 0.0);

      const double dpi_da = jacobi_p_derivative(a, i, 0.0, 0.0);
      const double dqj_db = jacobi_p_derivative(b, j, 2.0 * i + 1.0, 0.0);
      double dwarp_db = 0.0;
      if (i > 0)
        dwarp_db =
            -0.5 * static_cast<double>(i) * std::pow(half_one_minus_b, static_cast<double>(i - 1));

      const double dphi_da = norm * dpi_da * warp * qj;
      const double dphi_db = norm * (pi * dwarp_db * qj + pi * warp * dqj_db);

      gradients_xi[index] = dphi_da * da_dxi;
      gradients_eta[index] = dphi_da * da_deta + dphi_db;
      ++index;
    }
  }
}

} // namespace wasl
