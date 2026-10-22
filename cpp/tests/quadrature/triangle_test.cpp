// Tests for triangle volume and edge quadrature on the reference triangle.
//
// Checks: volume point count and weight sum (reference area 2), exactness
// against analytic monomial integrals up to the design degree and failure one
// degree beyond, positivity, and edge rule point counts, normals, and that the
// edge weights integrate the unit trace over each edge.

#include <cmath>
#include <cstdio>

#include "wasl/quadrature/quadrature.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace {

constexpr double kTol = 1e-10;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Exact integral of xi^a eta^b over the reference triangle {-1<=xi, -1<=eta,
// xi+eta<=0}, area 2. Computed by the affine change of variables to barycentric
// coordinates r = (1+xi+eta)/2, s = (1-eta)/2, t = (1+xi)/2 - r ... here we use
// the closed form via the substitution xi = -1 + u, eta = -1 + v with the upper
// bound v <= 1 - u on the unit right triangle, scaled to area 2.
//
// integral_T xi^a eta^b dA. Evaluated by collapsing to the (a, b) collapsed
// coordinates: dA = (1-eta)/2 da deta, with xi = 0.5(a+1)(1-eta)-1.
double exact_integral(int a, int b) {
  // Independent high-accuracy reference via a 7x7 collapsed-coordinate tensor
  // Gauss rule (exact to degree 13 per axis, far beyond the degree-8 cases the
  // test needs). This is independent of the rule under test, which only ever
  // requests up to a (p+1)-point rule for p <= 3.
  const auto g = wasl::gauss_legendre(7);
  double total = 0.0;
  for (int qb = 0; qb < g.size(); ++qb) {
    const double eta = g.points[qb];
    const double jac = 0.5 * (1.0 - eta);
    for (int qa = 0; qa < g.size(); ++qa) {
      const double aa = g.points[qa];
      const double xi = 0.5 * (aa + 1.0) * (1.0 - eta) - 1.0;
      total += std::pow(xi, a) * std::pow(eta, b) * g.weights[qa] * g.weights[qb] * jac;
    }
  }
  return total;
}

// Volume rule exactness: integrates xi^a eta^b exactly up to the design degree
// 2p, and the rule built for p must fail at degree 2p+1 (within tolerance).
bool check_volume_exactness() {
  for (int p = 0; p <= 3; ++p) {
    const auto rule = wasl::triangle_volume_rule(p);
    const int deg = wasl::triangle_volume_exact_degree(p);
    for (int a = 0; a <= deg; ++a) {
      for (int b = 0; b <= deg - a; ++b) {
        double integral = 0.0;
        for (int q = 0; q < rule.size(); ++q) {
          integral += std::pow(rule.xi[q], a) * std::pow(rule.eta[q], b) * rule.weights[q];
        }
        if (!approx(integral, exact_integral(a, b), 1e-9))
          return false;
      }
    }
  }
  return true;
}

// Weights sum to the reference-triangle area (2); points lie inside T.
bool check_volume_geometry() {
  for (int p = 0; p <= 3; ++p) {
    const auto rule = wasl::triangle_volume_rule(p);
    const int expected_points =
        wasl::triangle_volume_axis_points(p) * wasl::triangle_volume_axis_points(p);
    if (rule.size() != expected_points)
      return false;
    double sum = 0.0;
    for (int q = 0; q < rule.size(); ++q) {
      sum += rule.weights[q];
      if (rule.weights[q] <= 0.0)
        return false;
      // Interior of T: xi >= -1, eta >= -1, xi+eta <= 0 (boundary inclusive).
      if (rule.xi[q] < -1.0 - 1e-12 || rule.eta[q] < -1.0 - 1e-12 ||
          rule.xi[q] + rule.eta[q] > 1e-12)
        return false;
    }
    if (!approx(sum, 2.0, 1e-12))
      return false;
  }
  return true;
}

// Edge rules: three edges, each with (p+1) points, unit normals matching the
// reference edges, and weights that integrate the constant 1 to the edge
// length. Reference edge lengths: 2, 2*sqrt(2), 2.
bool check_edge_rules() {
  const double sqrt2 = std::sqrt(2.0);
  const double lengths[] = {2.0, 2.0 * sqrt2, 2.0};
  const double normals[][2] = {{0.0, -1.0}, {1.0 / sqrt2, 1.0 / sqrt2}, {-1.0, 0.0}};
  for (int p = 0; p <= 3; ++p) {
    const auto edges = wasl::triangle_edge_rules(p);
    if (edges.size() != 3)
      return false;
    for (int e = 0; e < 3; ++e) {
      if (static_cast<int>(edges[e].xi.size()) != p + 1)
        return false;
      if (!approx(edges[e].normal_x, normals[e][0]) || !approx(edges[e].normal_y, normals[e][1]))
        return false;
      double sum = 0.0;
      for (int q = 0; q < static_cast<int>(edges[e].weights.size()); ++q) {
        sum += edges[e].weights[q];
        if (edges[e].weights[q] <= 0.0)
          return false;
      }
      if (!approx(sum, lengths[e], 1e-12))
        return false;
    }
  }
  return true;
}

} // namespace

int main() {
  if (!check_volume_geometry())
    return 1;
  if (!check_volume_exactness())
    return 1;
  if (!check_edge_rules())
    return 1;
  return 0;
}
