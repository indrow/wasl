// Quadrature rules on the reference line xi in [-1, 1].
//
// Volume integrals use Gauss-Legendre rules. An n-point rule integrates
// polynomials up to degree 2 n - 1 exactly. Rules for n = 1..7 are tabulated
// to full double precision, covering every polynomial order Wasl targets in 1D.
//
// Endpoint (face) quadrature in 1D is degenerate: a face is a single point,
// so its weight is one. The two faces of a line element sit at xi = -1 (left,
// outward normal -1) and xi = +1 (right, outward normal +1). The endpoint rule
// is exposed so downstream code shares one face-evaluation shape with the
// later 2D edge rules.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace wasl {

struct QuadratureRule {
  std::vector<double> points;
  std::vector<double> weights;

  int size() const noexcept { return static_cast<int>(points.size()); }
};

namespace detail {

// One symmetric half of each Gauss-Legendre rule: the nonnegative node and its
// weight. A rule with a node at zero stores it once with weight w0. A rule with
// no node at zero stores node/weight pairs for the positive nodes; the mirror
// node shares the same weight.
struct SymmetricHalf {
  double node;
  double weight;
};

// Rules indexed by point count n (1..7). rules[n].node == 0.0 marks a rule that
// includes the origin; that entry is emitted once.
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
  if (n < 1 || n >= static_cast<int>(table.size())) {
    throw std::out_of_range("gauss_legendre: unsupported point count");
  }
  return table[n];
}

} // namespace detail

// Build the n-point Gauss-Legendre rule on [-1, 1].
inline QuadratureRule gauss_legendre(int n) {
  const auto &halves = detail::gauss_legendre_halves(n);
  QuadratureRule rule;
  for (auto it = halves.rbegin(); it != halves.rend(); ++it) {
    if (it->node == 0.0)
      continue; // origin emitted in the middle below
    rule.points.push_back(-it->node);
    rule.weights.push_back(it->weight);
  }
  for (const auto &h : halves) {
    rule.points.push_back(h.node);
    rule.weights.push_back(h.weight);
  }
  return rule;
}

// Highest polynomial degree integrated exactly by an n-point Gauss rule.
inline constexpr int gauss_legendre_exact_degree(int n) noexcept { return 2 * n - 1; }

// Endpoint face quadrature: the two faces of a reference line element.
// points are xi = -1 (left face) and xi = +1 (right face), each with weight one.
inline QuadratureRule endpoint_faces() {
  return QuadratureRule{
      std::vector<double>{-1.0, 1.0},
      std::vector<double>{1.0, 1.0},
  };
}

} // namespace wasl
