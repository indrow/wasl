// Tests for the orthonormal Dubiner basis on the reference triangle.
//
// Checks: mode count (P0-P3), constant-mode exactness, vertex values against
// the known mode-zero constant, discrete orthonormality through P3 via the
// triangle volume rule, and gradient consistency against centered finite
// differences.

#include <cmath>
#include <cstdio>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace {

constexpr double kTol = 1e-10;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Mode count matches triangle_mode_count at P0-P3.
bool check_mode_count() {
  for (int p = 0; p <= 3; ++p) {
    const std::size_t expected = wasl::triangle_mode_count(p);
    double vals[10];
    wasl::dubiner_basis(-0.5, -0.5, p, vals);
    if (expected != (p + 1) * (p + 2) / 2)
      return false;
    (void)vals;
  }
  return true;
}

// Mode zero is the constant 1/sqrt(2): the only mode with nonzero mean.
bool check_constant_mode() {
  const double pts[][2] = {{-1.0, -1.0}, {0.0, -1.0}, {-1.0, 0.0}, {-0.3, -0.4}, {-0.1, -0.1}};
  double vals[10];
  for (const auto &pt : pts) {
    wasl::dubiner_basis(pt[0], pt[1], 3, vals);
    if (!approx(vals[0], 1.0 / std::sqrt(2.0)))
      return false;
  }
  return true;
}

// Discrete orthonormality through P3: integrate phi_ij phi_kl over the
// reference triangle with the degree-8 volume rule (exact far beyond the
// degree-6 P3*P3 integrand). The mass matrix must be the identity.
bool check_orthonormality() {
  constexpr int p = 3;
  const std::size_t n = wasl::triangle_mode_count(p);
  const auto rule = wasl::triangle_volume_rule(4); // degree 8, well above 2*3

  std::vector<std::vector<double>> mass(n, std::vector<double>(n, 0.0));
  std::vector<double> vi(n);
  std::vector<double> vj(n);
  for (int q = 0; q < rule.size(); ++q) {
    wasl::dubiner_basis(rule.xi[q], rule.eta[q], p, vi.data());
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i; j < n; ++j) { // symmetric
        mass[i][j] += vi[i] * vi[j] * rule.weights[q];
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!approx(mass[i][i], 1.0))
      return false;
    for (std::size_t j = i + 1; j < n; ++j) {
      if (!approx(mass[i][j], 0.0))
        return false;
      if (!approx(mass[j][i], 0.0))
        return false;
    }
  }
  return true;
}

// The weights of the volume rule sum to the reference-triangle area (2).
bool check_volume_rule_area() {
  for (int p = 0; p <= 3; ++p) {
    const auto rule = wasl::triangle_volume_rule(p);
    double sum = 0.0;
    for (double w : rule.weights)
      sum += w;
    if (!approx(sum, 2.0, 1e-12))
      return false;
  }
  return true;
}

// Gradient of the constant mode is zero; gradients match centered finite
// differences at several interior points.
bool check_gradient_finite_difference() {
  constexpr double h = 1e-6;
  constexpr int p = 3;
  const std::size_t n = wasl::triangle_mode_count(p);
  const double pts[][2] = {{-0.8, -0.1}, {-0.3, -0.5}, {-0.2, -0.2}, {-0.6, -0.3}};
  std::vector<double> gx(n), gy(n), vp(n), vm(n);
  for (const auto &pt : pts) {
    wasl::dubiner_basis_gradient(pt[0], pt[1], p, gx.data(), gy.data());

    // Constant-mode gradient must be exactly zero.
    if (!approx(gx[0], 0.0, 1e-12) || !approx(gy[0], 0.0, 1e-12))
      return false;

    wasl::dubiner_basis(pt[0] + h, pt[1], p, vp.data());
    wasl::dubiner_basis(pt[0] - h, pt[1], p, vm.data());
    for (std::size_t m = 0; m < n; ++m) {
      const double fd = (vp[m] - vm[m]) / (2.0 * h);
      if (!approx(gx[m], fd, 1e-4))
        return false;
    }

    wasl::dubiner_basis(pt[0], pt[1] + h, p, vp.data());
    wasl::dubiner_basis(pt[0], pt[1] - h, p, vm.data());
    for (std::size_t m = 0; m < n; ++m) {
      const double fd = (vp[m] - vm[m]) / (2.0 * h);
      if (!approx(gy[m], fd, 1e-4))
        return false;
    }
  }
  return true;
}

// Reproduction of a constant field: projecting f(x,y) = c onto the basis gives
// only the mode-zero coefficient c * sqrt(2) (since phi_0 = 1/sqrt(2) and the
// cell average is c). All other modal moments must vanish.
bool check_constant_projection() {
  const auto rule = wasl::triangle_volume_rule(3);
  constexpr int p = 2;
  const std::size_t n = wasl::triangle_mode_count(p);
  const double c = 2.5;
  std::vector<double> moments(n, 0.0);
  std::vector<double> vals(n);
  for (int q = 0; q < rule.size(); ++q) {
    wasl::dubiner_basis(rule.xi[q], rule.eta[q], p, vals.data());
    for (std::size_t m = 0; m < n; ++m)
      moments[m] += c * vals[m] * rule.weights[q];
  }
  if (!approx(moments[0], c * std::sqrt(2.0)))
    return false;
  for (std::size_t m = 1; m < n; ++m) {
    if (!approx(moments[m], 0.0))
      return false;
  }
  return true;
}

} // namespace

int main() {
  if (!check_mode_count())
    return 1;
  if (!check_constant_mode())
    return 1;
  if (!check_orthonormality())
    return 1;
  if (!check_volume_rule_area())
    return 1;
  if (!check_gradient_finite_difference())
    return 1;
  if (!check_constant_projection())
    return 1;
  return 0;
}
