// Tests for the orthonormal Legendre basis on the reference line.
//
// Checks: mode count (P0-P3), constant-mode exactness, endpoint values against
// known closed forms, and discrete orthonormality through P3.

#include <cmath>

#include "wasl/basis/legendre.hpp"
#include "wasl/quadrature/quadrature.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Exact endpoint values of orthonormal Legendre polynomials phi_j(1) and
// phi_j(-1) = (-1)^j phi_j(1), with phi_j(1) = sqrt((2j+1)/2).
bool check_endpoints() {
  // phi_j(+1) = sqrt((2j+1)/2); phi_j(-1) = (-1)^j phi_j(+1).
  double vals[4];
  for (int p = 0; p <= 3; ++p) {
    wasl::legendre_basis(1.0, p, vals);
    for (int j = 0; j <= p; ++j) {
      if (!approx(vals[j], wasl::legendre_norm(j)))
        return false;
    }
    wasl::legendre_basis(-1.0, p, vals);
    for (int j = 0; j <= p; ++j) {
      const double expected = ((j % 2 == 0) ? 1.0 : -1.0) * wasl::legendre_norm(j);
      if (!approx(vals[j], expected))
        return false;
    }
  }
  return true;
}

// Mode zero is constant 1/sqrt(2) at several points.
bool check_constant_mode() {
  double vals[4];
  const double xs[] = {-1.0, -0.4, 0.0, 0.37, 1.0};
  for (double xi : xs) {
    wasl::legendre_basis(xi, 3, vals);
    if (!approx(vals[0], 1.0 / std::sqrt(2.0)))
      return false;
  }
  return true;
}

// Mode one is the linear mode sqrt(3/2) * xi.
bool check_linear_mode() {
  double vals[2];
  const double xs[] = {-1.0, -0.3, 0.0, 0.5, 1.0};
  for (double xi : xs) {
    wasl::legendre_basis(xi, 1, vals);
    if (!approx(vals[1], std::sqrt(1.5) * xi))
      return false;
  }
  return true;
}

// Discrete orthonormality: integrate phi_i phi_j with a high-order Gauss rule.
// The mass matrix M_ij must be the identity.
bool check_orthonormality() {
  // 7-point rule integrates polynomials up to degree 13 exactly; P3*P3 is
  // degree 6, so this is far more than enough.
  const auto rule = wasl::gauss_legendre(7);
  constexpr int p = 3;
  constexpr int n = p + 1;
  double mass[n][n] = {{0}};
  double vi[n];
  for (int q = 0; q < rule.size(); ++q) {
    wasl::legendre_basis(rule.points[q], p, vi);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        mass[i][j] += vi[i] * vi[j] * rule.weights[q];
  }
  for (int i = 0; i < n; ++i) {
    if (!approx(mass[i][i], 1.0, 1e-10))
      return false;
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;
      if (!approx(mass[i][j], 0.0, 1e-10))
        return false;
    }
  }
  return true;
}

// Gradient of the constant mode is zero; gradient of mode 1 is sqrt(3/2).
bool check_gradient_basic() {
  double grads[4];
  wasl::legendre_basis_gradient(0.42, 3, grads);
  if (!approx(grads[0], 0.0))
    return false;
  if (!approx(grads[1], std::sqrt(1.5)))
    return false;
  return true;
}

// Gradient consistency: numerical finite difference of the basis matches the
// analytic gradient.
bool check_gradient_finite_difference() {
  constexpr double h = 1e-6;
  constexpr int p = 3;
  double vp[p + 1], vm[p + 1], grads[p + 1];
  const double xs[] = {-0.9, -0.1, 0.25, 0.8};
  for (double xi : xs) {
    wasl::legendre_basis(xi + h, p, vp);
    wasl::legendre_basis(xi - h, p, vm);
    wasl::legendre_basis_gradient(xi, p, grads);
    for (int j = 0; j <= p; ++j) {
      const double fd = (vp[j] - vm[j]) / (2.0 * h);
      if (!approx(grads[j], fd, 1e-4))
        return false;
    }
  }
  return true;
}

} // namespace

int main() {
  if (!check_endpoints())
    return 1;
  if (!check_constant_mode())
    return 1;
  if (!check_linear_mode())
    return 1;
  if (!check_orthonormality())
    return 1;
  if (!check_gradient_basic())
    return 1;
  if (!check_gradient_finite_difference())
    return 1;
  return 0;
}
