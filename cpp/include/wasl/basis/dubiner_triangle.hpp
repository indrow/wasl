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
// where P_n^(alpha,beta) is the standard (orthogonal) Jacobi polynomial. With
// this normalization
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

#include "wasl/basis/reference_element.hpp"

namespace wasl {
namespace detail {

// Standard Jacobi polynomial P_n^(alpha, beta)(x), orthogonal with respect to
// the weight (1 - x)^alpha (1 + x)^beta on [-1, 1]. Base cases P_0 = 1 and
// P_1 = ((alpha - beta) + (alpha + beta + 2) x) / 2, then the three-term
// recurrence (Abramowitz-Stegner 22.7.1). The recurrence coefficients stay
// nonzero for the nonnegative alpha, beta used here.
inline double jacobi_p(double x, int n, double alpha, double beta) noexcept {
  if (n <= 0)
    return 1.0;
  double p_prev = 1.0;                                             // P_0
  double p_curr = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x); // P_1
  for (int k = 1; k < n; ++k) {
    const double ab = 2.0 * k + alpha + beta;
    const double a = 2.0 * (k + 1) * (k + alpha + beta + 1.0) * ab;
    const double b = (ab + 1.0) * (ab * (ab + 2.0) * x + (alpha * alpha - beta * beta));
    const double c = 2.0 * (k + alpha) * (k + beta) * (ab + 2.0);
    const double p_next = (b * p_curr - c * p_prev) / a;
    p_prev = p_curr;
    p_curr = p_next;
  }
  return p_curr;
}

// Derivative d/dx of the standard Jacobi polynomial, via the identity
//   d/dx P_n^(alpha,beta)(x) = (n + alpha + beta + 1) / 2 * P_{n-1}^(alpha+1,beta+1)(x).
// Returns 0 for the degree-zero mode.
inline double jacobi_p_derivative(double x, int n, double alpha, double beta) noexcept {
  if (n <= 0)
    return 0.0;
  return 0.5 * (n + alpha + beta + 1.0) * jacobi_p(x, n - 1, alpha + 1.0, beta + 1.0);
}

} // namespace detail

// Evaluate orthonormal Dubiner modes 0..(triangle_mode_count(p)-1) at (xi, eta)
// into values. values must point to at least triangle_mode_count(p) doubles.
inline void dubiner_basis(double xi, double eta, int p, double *values) noexcept {
  const double b = eta;
  const double a = (eta < 1.0) ? (2.0 + 2.0 * xi) / (1.0 - eta) - 1.0 : -1.0;
  const double half_one_minus_b = 0.5 * (1.0 - b);

  int idx = 0;
  for (int n = 0; n <= p; ++n) {
    for (int i = 0; i <= n; ++i) {
      const int j = n - i;
      const double norm = std::sqrt((2.0 * i + 1.0) * (n + 1.0) / 2.0);
      const double pi = detail::jacobi_p(a, i, 0.0, 0.0);
      const double warp = std::pow(half_one_minus_b, static_cast<double>(i));
      const double qj = detail::jacobi_p(b, j, 2.0 * i + 1.0, 0.0);
      values[idx++] = norm * pi * warp * qj;
    }
  }
}

// Evaluate orthonormal Dubiner gradients (d/dxi, d/deta) of modes 0..count-1 at
// (xi, eta) into grads_xi and grads_eta. Each must point to at least
// triangle_mode_count(p) doubles. Gradients are defined for interior points
// (eta < 1); the apex eta = 1 is a collapsed-coordinate singularity.
inline void dubiner_basis_gradient(double xi, double eta, int p, double *grads_xi,
                                   double *grads_eta) noexcept {
  const double b = eta;
  const double a = (eta < 1.0) ? (2.0 + 2.0 * xi) / (1.0 - eta) - 1.0 : -1.0;
  const double half_one_minus_b = 0.5 * (1.0 - b);
  const double inv_one_minus_b = 1.0 / (1.0 - b); // finite for eta < 1

  // Collapsed-coordinate chain-rule factors.
  const double da_dxi = 2.0 * inv_one_minus_b;
  const double da_deta = (a + 1.0) * inv_one_minus_b;

  int idx = 0;
  for (int n = 0; n <= p; ++n) {
    for (int i = 0; i <= n; ++i) {
      const int j = n - i;
      const double norm = std::sqrt((2.0 * i + 1.0) * (n + 1.0) / 2.0);
      const double pi = detail::jacobi_p(a, i, 0.0, 0.0);
      const double warp = std::pow(half_one_minus_b, static_cast<double>(i));
      const double qj = detail::jacobi_p(b, j, 2.0 * i + 1.0, 0.0);

      const double dpi_da = detail::jacobi_p_derivative(a, i, 0.0, 0.0);
      const double dqj_db = detail::jacobi_p_derivative(b, j, 2.0 * i + 1.0, 0.0);
      double dwarp_db = 0.0;
      if (i > 0)
        dwarp_db =
            -0.5 * static_cast<double>(i) * std::pow(half_one_minus_b, static_cast<double>(i - 1));

      // d phi / d a and d phi / d b in collapsed coordinates.
      const double dphi_da = norm * dpi_da * warp * qj;
      const double dphi_db = norm * (pi * dwarp_db * qj + pi * warp * dqj_db);

      grads_xi[idx] = dphi_da * da_dxi;
      grads_eta[idx] = dphi_da * da_deta + dphi_db;
      ++idx;
    }
  }
}

} // namespace wasl
