// Volume quadrature on the reference triangle
//   {(xi, eta) : xi >= -1, eta >= -1, xi + eta <= 0},  area = 2.
//
// Rules are the collapsed-coordinate tensor product of two Gauss rules, the
// same coordinate map the Dubiner basis uses:
//   b = eta,  a = (2 + 2 xi) / (1 - eta) - 1,
// with Jacobian da db -> dxi deta of (1 - eta) / 2. Pairing two (p+1)-point
// rules integrates the degree-2p mass and stiffness integrands of a Pp
// triangle exactly.

#pragma once

#include <vector>

#include "wasl/quadrature/gauss_legendre.hpp"

namespace wasl {

struct TriangleQuadratureRule {
  std::vector<double> xi;  // [point] reference xi
  std::vector<double> eta; // [point] reference eta
  std::vector<double> weights;

  int size() const noexcept { return static_cast<int>(xi.size()); }
};

inline constexpr int triangle_volume_exact_degree(int p) noexcept { return 2 * p; }

// Gauss points along each collapsed axis of the order-p rule.
inline int triangle_volume_axis_points(int p) noexcept { return p + 1; }

// The weights include the collapsed-coordinate Jacobian, so
// sum_q w_q f(xi_q, eta_q) approximates integral_T f dA.
inline TriangleQuadratureRule triangle_volume_rule(int p) {
  const int n = triangle_volume_axis_points(p);
  const QuadratureRule b_rule = gauss_legendre(n);
  const QuadratureRule a_rule = gauss_legendre(n);

  TriangleQuadratureRule rule;
  for (int qb = 0; qb < b_rule.size(); ++qb) {
    const double eta = b_rule.points[qb];
    const double jacobian = 0.5 * (1.0 - eta);
    for (int qa = 0; qa < a_rule.size(); ++qa) {
      const double a = a_rule.points[qa];
      rule.xi.push_back(0.5 * (a + 1.0) * (1.0 - eta) - 1.0);
      rule.eta.push_back(eta);
      rule.weights.push_back(b_rule.weights[qb] * a_rule.weights[qa] * jacobian);
    }
  }
  return rule;
}

} // namespace wasl
