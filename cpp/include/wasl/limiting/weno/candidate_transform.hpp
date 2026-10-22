#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/mesh/triangle_geometry.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace wasl {

// L2-project a polynomial from its affine source triangle into the target
// triangle basis. Both coefficient arrays contain triangle_mode_count(p)
// entries in the documented Dubiner ordering.
inline void project_neighbor_polynomial(const TriangleMap &source, const TriangleMap &target,
                                        int polynomial_order, const double *source_coefficients,
                                        double *target_coefficients) {
  const std::size_t mode_count = triangle_mode_count(polynomial_order);
  const TriangleQuadratureRule quadrature = triangle_volume_rule(polynomial_order);
  std::vector<double> source_basis(mode_count);
  std::vector<double> target_basis(mode_count);
  for (std::size_t mode = 0; mode < mode_count; ++mode)
    target_coefficients[mode] = 0.0;

  for (int q = 0; q < quadrature.size(); ++q) {
    double x = 0.0, y = 0.0;
    target.reference_to_physical(quadrature.xi[q], quadrature.eta[q], x, y);
    double source_xi = 0.0, source_eta = 0.0;
    source.physical_to_reference(x, y, source_xi, source_eta);
    dubiner_basis(source_xi, source_eta, polynomial_order, source_basis.data());
    dubiner_basis(quadrature.xi[q], quadrature.eta[q], polynomial_order, target_basis.data());

    double value = 0.0;
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      value += source_coefficients[mode] * source_basis[mode];
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      target_coefficients[mode] += quadrature.weights[q] * value * target_basis[mode];
  }
}

inline void make_weno_candidate(const TriangleMap &source, const TriangleMap &target,
                                int polynomial_order, const double *source_coefficients,
                                double target_mean_mode, double *candidate) {
  project_neighbor_polynomial(source, target, polynomial_order, source_coefficients, candidate);
  candidate[0] = target_mean_mode;
}

} // namespace wasl
