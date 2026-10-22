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
// values must point to at least p + 1 doubles.
inline void legendre_basis(double xi, int p, double *values) noexcept {
  values[0] = legendre_norm(0); // phi_0 = 1 / sqrt(2), constant
  if (p <= 0)
    return;

  // Standard Legendre recurrence P_{j+1} = ((2j+1) xi P_j - j P_{j-1}) / (j+1).
  double pjm1 = 1.0; // P_{j-1}
  double pj = xi;    // P_j
  values[1] = legendre_norm(1) * pj;
  for (int j = 1; j < p; ++j) {
    const double pj_next = ((2.0 * j + 1.0) * xi * pj - j * pjm1) / (j + 1);
    pjm1 = pj;
    pj = pj_next;
    values[j + 1] = legendre_norm(j + 1) * pj;
  }
}

// Evaluate orthonormal gradients d/dxi of modes 0..p at xi into grads[0..p].
// grads must point to at least p + 1 doubles.
inline void legendre_basis_gradient(double xi, int p, double *grads) noexcept {
  grads[0] = 0.0; // derivative of the constant mode
  if (p <= 0)
    return;

  // Coupled recurrence on the derivatives of the standard Legendre polynomials:
  //   (j+1) P'_{j+1} = (2j+1)(P_j + xi P'_j) - j P'_{j-1}.
  double pjm1 = 1.0, pj = xi;    // P_{j-1}, P_j
  double dpjm1 = 0.0, dpj = 1.0; // P'_{j-1}, P'_j
  grads[1] = legendre_norm(1) * dpj;
  for (int j = 1; j < p; ++j) {
    const double pj_next = ((2.0 * j + 1.0) * xi * pj - j * pjm1) / (j + 1);
    const double dpj_next = ((2.0 * j + 1.0) * (pj + xi * dpj) - j * dpjm1) / (j + 1);
    pjm1 = pj;
    pj = pj_next;
    dpjm1 = dpj;
    dpj = dpj_next;
    grads[j + 1] = legendre_norm(j + 1) * dpj;
  }
}

// Combined values and gradients at xi. Either pointer may be null.
inline void legendre_basis_with_gradient(double xi, int p, double *values, double *grads) noexcept {
  if (values != nullptr)
    legendre_basis(xi, p, values);
  if (grads != nullptr)
    legendre_basis_gradient(xi, p, grads);
}

} // namespace wasl
