// Gauss-Legendre rules on the reference line xi in [-1, 1].
//
// An n-point rule integrates polynomials up to degree 2 n - 1 exactly. Rules for
// n = 1..7 are tabulated to full double precision, covering every polynomial
// order Wasl targets in 1D.

#pragma once

#include <stdexcept>
#include <vector>

#include "wasl/quadrature/quadrature_rule.hpp"

namespace wasl {
namespace detail {

// The nonnegative half of a symmetric rule. A node at the origin is stored
// once; every positive node implies its mirror with the same weight.
struct SymmetricHalf {
  double node;
  double weight;
};

inline const std::vector<SymmetricHalf> &gauss_legendre_halves(int n) {
  static const std::vector<std::vector<SymmetricHalf>> table = {
      {},                                                                    // n = 0 unused
      {{0.0, 2.0}},                                                          // n = 1
      {{0.5773502691896257, 1.0}},                                           // n = 2
      {{0.0, 0.8888888888888888}, {0.7745966692414834, 0.5555555555555556}}, // n = 3
      {{0.33998104358485626, 0.6521451548625461},
       {0.8611363115940526, 0.34785484513745385}}, // n = 4
      {{0.0, 0.5688888888888889},
       {0.5384693101056831, 0.47862867049936647},
       {0.9061798459386640, 0.2369268850561891}}, // n = 5
      {{0.2386191860831969, 0.4679139345726910},
       {0.6612093864662645, 0.3607615730481386},
       {0.9324695142031521, 0.17132449237917036}}, // n = 6
      {{0.0, 0.4179591836734694},
       {0.4058451513773972, 0.3818300505051188},
       {0.7415311855993945, 0.2797053914892766},
       {0.9491079123427585, 0.1294849661688697}}, // n = 7
  };
  if (n < 1 || n >= static_cast<int>(table.size()))
    throw std::out_of_range("gauss_legendre: unsupported point count");
  return table[n];
}

} // namespace detail

// The n-point Gauss-Legendre rule on [-1, 1], with points ascending.
inline QuadratureRule gauss_legendre(int n) {
  const auto &halves = detail::gauss_legendre_halves(n);
  QuadratureRule rule;
  for (auto half = halves.rbegin(); half != halves.rend(); ++half) {
    if (half->node == 0.0)
      continue; // the origin is emitted once, below
    rule.points.push_back(-half->node);
    rule.weights.push_back(half->weight);
  }
  for (const auto &half : halves) {
    rule.points.push_back(half.node);
    rule.weights.push_back(half.weight);
  }
  return rule;
}

inline constexpr int gauss_legendre_exact_degree(int n) noexcept { return 2 * n - 1; }

} // namespace wasl
