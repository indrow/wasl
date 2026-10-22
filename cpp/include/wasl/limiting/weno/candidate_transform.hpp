// Neighbor candidate polynomials for WENO on triangles.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/modal_expansion.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/mesh/triangle/triangle_geometry.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {

// L2-project the polynomial of an affine source triangle into the target
// triangle basis. Both arrays hold triangle_mode_count(p) Dubiner coefficients.
inline void project_neighbor_polynomial(const TriangleMap &source, const TriangleMap &target,
                                        int polynomial_order, const double *source_coefficients,
                                        double *target_coefficients) {
  const std::size_t mode_count = triangle_mode_count(polynomial_order);
  const TriangleQuadratureRule quadrature = triangle_volume_rule(polynomial_order);
  std::vector<double> source_basis(mode_count), target_basis(mode_count);
  for (std::size_t mode = 0; mode < mode_count; ++mode)
    target_coefficients[mode] = 0.0;

  for (int q = 0; q < quadrature.size(); ++q) {
    double x = 0.0, y = 0.0;
    target.reference_to_physical(quadrature.xi[q], quadrature.eta[q], x, y);
    double source_xi = 0.0, source_eta = 0.0;
    source.physical_to_reference(x, y, source_xi, source_eta);
    dubiner_basis(source_xi, source_eta, polynomial_order, source_basis.data());
    dubiner_basis(quadrature.xi[q], quadrature.eta[q], polynomial_order, target_basis.data());

    const double value = modal_sum(source_coefficients, source_basis.data(), mode_count);
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      target_coefficients[mode] += quadrature.weights[q] * value * target_basis[mode];
  }
}

// A projected neighbor with the target's own mean, so blending preserves the
// cell average.
inline void make_weno_candidate(const TriangleMap &source, const TriangleMap &target,
                                int polynomial_order, const double *source_coefficients,
                                double target_mean_mode, double *candidate) {
  project_neighbor_polynomial(source, target, polynomial_order, source_coefficients, candidate);
  candidate[0] = target_mean_mode;
}

} // namespace wasl
