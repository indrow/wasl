// Tests for Gauss-Legendre and endpoint quadrature on the reference line.
//
// Checks: point count, symmetry, weights sum to 2, exact integration of
// monomials up to the rule's exactness degree, and the endpoint rule shape.

#include <cmath>

#include "wasl/quadrature/quadrature.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Weights of a Gauss-Legendre rule sum to 2 on [-1, 1] (integral of 1).
bool check_weight_sum() {
  for (int n = 1; n <= 7; ++n) {
    const auto rule = wasl::gauss_legendre(n);
    double sum = 0.0;
    for (double w : rule.weights)
      sum += w;
    if (!approx(sum, 2.0))
      return false;
  }
  return true;
}

// Every rule must have exactly n points.
bool check_point_count() {
  for (int n = 1; n <= 7; ++n) {
    const auto rule = wasl::gauss_legendre(n);
    if (rule.size() != n)
      return false;
  }
  return true;
}

// xi^k integrates exactly for k = 0..(2n-1) and fails at k = 2n.
bool check_exactness() {
  for (int n = 1; n <= 7; ++n) {
    const auto rule = wasl::gauss_legendre(n);
    const int exact_deg = wasl::gauss_legendre_exact_degree(n);
    for (int k = 0; k <= exact_deg; ++k) {
      double integral = 0.0;
      for (int q = 0; q < rule.size(); ++q) {
        integral += std::pow(rule.points[q], k) * rule.weights[q];
      }
      // integral of xi^k over [-1,1] = 0 for odd k, 2/(k+1) for even k.
      const double exact = (k % 2 == 1) ? 0.0 : 2.0 / (k + 1.0);
      if (!approx(integral, exact, 1e-11))
        return false;
    }
    // One degree beyond exactness must NOT integrate exactly.
    const int k = exact_deg + 1;
    double integral = 0.0;
    for (int q = 0; q < rule.size(); ++q) {
      integral += std::pow(rule.points[q], k) * rule.weights[q];
    }
    const double exact = (k % 2 == 1) ? 0.0 : 2.0 / (k + 1.0);
    if (approx(integral, exact, 1e-6))
      return false;
  }
  return true;
}

// Symmetry: points come in +/- pairs (ignoring the origin), equal weights.
bool check_symmetry() {
  for (int n = 1; n <= 7; ++n) {
    const auto rule = wasl::gauss_legendre(n);
    for (int i = 0; i < rule.size(); ++i) {
      bool found = false;
      for (int j = 0; j < rule.size(); ++j) {
        if (approx(rule.points[j], -rule.points[i]) && approx(rule.weights[j], rule.weights[i])) {
          found = true;
          break;
        }
      }
      if (!found)
        return false;
    }
  }
  return true;
}

// Endpoint rule: two faces at xi = -1, +1, each weight one.
bool check_endpoint_rule() {
  const auto rule = wasl::endpoint_faces();
  if (rule.size() != 2)
    return false;
  if (!approx(rule.points[0], -1.0))
    return false;
  if (!approx(rule.points[1], 1.0))
    return false;
  if (!approx(rule.weights[0], 1.0))
    return false;
  if (!approx(rule.weights[1], 1.0))
    return false;
  return true;
}

} // namespace

int main() {
  if (!check_weight_sum())
    return 1;
  if (!check_point_count())
    return 1;
  if (!check_exactness())
    return 1;
  if (!check_symmetry())
    return 1;
  if (!check_endpoint_rule())
    return 1;
  return 0;
}
