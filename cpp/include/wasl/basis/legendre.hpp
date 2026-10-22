// Orthonormal Legendre modal basis on the reference line xi in [-1, 1].
//
// Mode j is the degree-j orthonormal Legendre polynomial
//   phi_j(xi) = sqrt((2 j + 1) / 2) * P_j(xi)
// where P_j is the standard Legendre polynomial. With this normalization
//   integral_{-1}^{1} phi_i(xi) phi_j(xi) d xi = delta_ij.
// Mode zero is phi_0 = 1 / sqrt(2): the only mode with nonzero mean, so it
// carries the cell average. A constant field is represented by a single mode.
//
// These evaluators are header-only and allocation-free. Values and gradients
// use the coupled Bonnet recurrence, which is stable at xi = +/- 1.

#pragma once

#include <cmath>

namespace wasl {

// Normalization factor sqrt((2 j + 1) / 2) for mode j.
inline double legendre_norm(int j) noexcept { return std::sqrt((2.0 * j + 1.0) / 2.0); }

// Evaluate orthonormal modes 0..p at xi into values[0..p].
inline void legendre_basis(double xi, int p, double *values) noexcept {
  values[0] = legendre_norm(0);
  if (p <= 0)
    return;

  // P_{j+1} = ((2j+1) xi P_j - j P_{j-1}) / (j+1).
  double p_previous = 1.0;
  double p_current = xi;
  values[1] = legendre_norm(1) * p_current;
  for (int j = 1; j < p; ++j) {
    const double p_next = ((2.0 * j + 1.0) * xi * p_current - j * p_previous) / (j + 1);
    p_previous = p_current;
    p_current = p_next;
    values[j + 1] = legendre_norm(j + 1) * p_current;
  }
}

// Evaluate orthonormal gradients d/dxi of modes 0..p at xi into gradients[0..p].
inline void legendre_basis_gradient(double xi, int p, double *gradients) noexcept {
  gradients[0] = 0.0;
  if (p <= 0)
    return;

  // (j+1) P'_{j+1} = (2j+1)(P_j + xi P'_j) - j P'_{j-1}.
  double p_previous = 1.0, p_current = xi;
  double dp_previous = 0.0, dp_current = 1.0;
  gradients[1] = legendre_norm(1) * dp_current;
  for (int j = 1; j < p; ++j) {
    const double p_next = ((2.0 * j + 1.0) * xi * p_current - j * p_previous) / (j + 1);
    const double dp_next =
        ((2.0 * j + 1.0) * (p_current + xi * dp_current) - j * dp_previous) / (j + 1);
    p_previous = p_current;
    p_current = p_next;
    dp_previous = dp_current;
    dp_current = dp_next;
    gradients[j + 1] = legendre_norm(j + 1) * dp_current;
  }
}

} // namespace wasl
