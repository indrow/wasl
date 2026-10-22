// Quadrature rules on the reference triangle
//   {(xi, eta) : xi >= -1, eta >= -1, xi + eta <= 0},  area = 2,
// with vertices v0 = (-1, -1), v1 = (1, -1), v2 = (-1, 1).
//
// Volume rules are built as the collapsed-coordinate tensor product of two
// Gauss rules, the same coordinate map the Dubiner basis uses:
//   b = eta,  a = (2 + 2 xi) / (1 - eta) - 1,
// with Jacobian da db -> dxi deta of (1 - eta) / 2. A degree-(2p) rule pairs an
// (p+1)-point rule in b with a (p+1)-point rule in a, integrating the degree-2p
// mass/stiffness integrands of a Pp triangle exactly.
//
// Edge rules map Gauss-Legendre points onto each reference edge and expose the
// reference outward normal, matching the 1D endpoint_faces shape so surface
// code shares one face-evaluation form with the line element.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/quadrature/quadrature.hpp"

namespace wasl {

struct TriangleQuadratureRule {
  std::vector<double> xi;  // [point] reference xi
  std::vector<double> eta; // [point] reference eta
  std::vector<double> weights;

  int size() const noexcept { return static_cast<int>(xi.size()); }
};

// Highest total polynomial degree integrated exactly by the volume rule built
// for polynomial order p. The collapsed tensor of two (p+1)-point Gauss rules
// integrates degree 2p exactly.
inline constexpr int triangle_volume_exact_degree(int p) noexcept { return 2 * p; }

// Number of (p+1)-point Gauss nodes used along each collapsed axis.
inline int triangle_volume_axis_points(int p) noexcept { return p + 1; }

// Build the degree-(2p) volume quadrature rule on the reference triangle. The
// weights already include the collapsed-coordinate Jacobian (1 - eta) / 2, so
// summing w_q f(xi_q, eta_q) approximates integral_T f dA.
inline TriangleQuadratureRule triangle_volume_rule(int p) {
  const int n = triangle_volume_axis_points(p);
  const QuadratureRule b_rule = gauss_legendre(n); // eta axis
  const QuadratureRule a_rule = gauss_legendre(n); // collapsed a axis

  TriangleQuadratureRule rule;
  for (int qb = 0; qb < b_rule.size(); ++qb) {
    const double eta = b_rule.points[qb];
    const double wb = b_rule.weights[qb];
    const double jac = 0.5 * (1.0 - eta); // |dxi deta / da db|
    for (int qa = 0; qa < a_rule.size(); ++qa) {
      const double a = a_rule.points[qa];
      const double wa = a_rule.weights[qa];
      const double xi = 0.5 * (a + 1.0) * (1.0 - eta) - 1.0;
      rule.xi.push_back(xi);
      rule.eta.push_back(eta);
      rule.weights.push_back(wb * wa * jac);
    }
  }
  return rule;
}

// A reference edge, parametrized by a 1D Gauss rule mapped onto the edge, with
// the reference outward unit normal. Edges follow the vertex pairs
//   edge 0: v0 -> v1,  the bottom edge eta = -1, outward normal (0, -1);
//   edge 1: v1 -> v2,  the hypotenuse xi + eta = 0, outward normal (1, 1)/sqrt2;
//   edge 2: v2 -> v0,  the left edge xi = -1, outward normal (-1, 0).
struct TriangleEdgeRule {
  std::vector<double> xi;      // [point] reference xi on the edge
  std::vector<double> eta;     // [point] reference eta on the edge
  std::vector<double> weights; // already scaled by the edge length
  double normal_x;             // reference outward unit normal x
  double normal_y;             // reference outward unit normal y
};

// Gauss-Legendre points on [-1, 1] scaled onto [x0, x1]: node t and weight w.
inline void map_segment(double x0, double x1, double t, double w, double &x, double &weight) {
  x = 0.5 * (x0 + x1) + 0.5 * (x1 - x0) * t;
  weight = 0.5 * (x1 - x0) * w;
}

// Build the edge quadrature rules for a Pp triangle. Each edge carries enough
// Gauss points to integrate the degree-(2p+1) surface trace product exactly.
inline std::vector<TriangleEdgeRule> triangle_edge_rules(int p) {
  const QuadratureRule g = gauss_legendre(p + 1);

  std::vector<TriangleEdgeRule> edges(3);
  // Edge 0: bottom, eta = -1, v0=(-1,-1) to v1=(1,-1).
  edges[0].normal_x = 0.0;
  edges[0].normal_y = -1.0;
  for (int q = 0; q < g.size(); ++q) {
    double x = 0.0, w = 0.0;
    map_segment(-1.0, 1.0, g.points[q], g.weights[q], x, w);
    edges[0].xi.push_back(x);
    edges[0].eta.push_back(-1.0);
    edges[0].weights.push_back(w);
  }

  // Edge 1: hypotenuse, xi + eta = 0, v1=(1,-1) to v2=(-1,1).
  // Length 2*sqrt(2); parametrize by reference coordinate s in [-1, 1].
  edges[1].normal_x = 1.0 / std::sqrt(2.0);
  edges[1].normal_y = 1.0 / std::sqrt(2.0);
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (int q = 0; q < g.size(); ++q) {
    double junk = 0.0, w = 0.0;
    map_segment(-1.0, 1.0, g.points[q], g.weights[q], junk, w);
    const double s = g.points[q];
    // v1 + 0.5*(s+1)*(v2-v1).
    const double xi = 1.0 + 0.5 * (s + 1.0) * (-2.0);
    const double eta = -1.0 + 0.5 * (s + 1.0) * (2.0);
    edges[1].xi.push_back(xi);
    edges[1].eta.push_back(eta);
    edges[1].weights.push_back(w * 2.0 * inv_sqrt2); // scale to physical length 2*sqrt2
  }

  // Edge 2: left, xi = -1, v2=(-1,1) to v0=(-1,-1). The quadrature nodes
  // parametrize eta in [-1, 1] with positive weights; the reference outward
  // normal carries the orientation, independent of node ordering.
  edges[2].normal_x = -1.0;
  edges[2].normal_y = 0.0;
  for (int q = 0; q < g.size(); ++q) {
    double y = 0.0, w = 0.0;
    map_segment(-1.0, 1.0, g.points[q], g.weights[q], y, w);
    edges[2].xi.push_back(-1.0);
    edges[2].eta.push_back(y);
    edges[2].weights.push_back(w);
  }
  return edges;
}

} // namespace wasl
