// Edge quadrature on the reference triangle with vertices
//   v0 = (-1, -1),  v1 = (1, -1),  v2 = (-1, 1).
//
// Edges follow the vertex pairs
//   edge 0: v0 -> v1,  the bottom edge eta = -1, outward normal (0, -1);
//   edge 1: v1 -> v2,  the hypotenuse xi + eta = 0, outward normal (1, 1)/sqrt2;
//   edge 2: v2 -> v0,  the left edge xi = -1, outward normal (-1, 0).
// Each edge carries p + 1 Gauss points, enough to integrate the degree-(2p+1)
// surface trace product of a Pp triangle exactly.

#pragma once

#include <cmath>
#include <vector>

#include "wasl/quadrature/gauss_legendre.hpp"

namespace wasl {

struct TriangleEdgeRule {
  std::vector<double> xi;      // [point] reference xi on the edge
  std::vector<double> eta;     // [point] reference eta on the edge
  std::vector<double> weights; // scaled by the reference edge length
  double normal_x;             // reference outward unit normal
  double normal_y;
};

inline std::vector<TriangleEdgeRule> triangle_edge_rules(int p) {
  const QuadratureRule gauss = gauss_legendre(p + 1);
  const double inverse_sqrt2 = 1.0 / std::sqrt(2.0);

  std::vector<TriangleEdgeRule> edges(3);
  edges[0].normal_x = 0.0;
  edges[0].normal_y = -1.0;
  edges[1].normal_x = inverse_sqrt2;
  edges[1].normal_y = inverse_sqrt2;
  edges[2].normal_x = -1.0;
  edges[2].normal_y = 0.0;
  for (int q = 0; q < gauss.size(); ++q) {
    const double s = gauss.points[q];
    const double weight = gauss.weights[q];

    edges[0].xi.push_back(s);
    edges[0].eta.push_back(-1.0);
    edges[0].weights.push_back(weight);

    // v1 + (s + 1) / 2 (v2 - v1), scaled to the hypotenuse length 2 sqrt2.
    edges[1].xi.push_back(1.0 + 0.5 * (s + 1.0) * (-2.0));
    edges[1].eta.push_back(-1.0 + 0.5 * (s + 1.0) * (2.0));
    edges[1].weights.push_back(weight * 2.0 * inverse_sqrt2);

    // Points ascend in eta; the outward normal carries the orientation.
    edges[2].xi.push_back(-1.0);
    edges[2].eta.push_back(s);
    edges[2].weights.push_back(weight);
  }
  return edges;
}

} // namespace wasl
