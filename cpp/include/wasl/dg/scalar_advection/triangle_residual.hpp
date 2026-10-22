// Modal DG residual for scalar advection on affine triangles.
//
// Velocity is constant or a spatial function. Coefficients use the scalar
// solver layout U[element][mode]. Boundary exterior traces use
// U[face][edge quadrature point]; interior entries are ignored.

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/modal_expansion.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/triangle_face_trace.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"
#include "wasl/quadrature/triangle_edge_rules.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {

// VelocityFunction: (x, y) -> ScalarAdvection2D.
template <typename VelocityFunction>
inline void triangle_scalar_advection_residual(const TriangleMesh &mesh, VelocityFunction velocity,
                                               int polynomial_order, const double *coefficients,
                                               const double *boundary_exterior_state, double *rhs) {
  const std::size_t mode_count = triangle_mode_count(polynomial_order);
  const int edge_point_count = polynomial_order + 1;
  const auto volume_rule = triangle_volume_rule(polynomial_order);
  const auto edge_rules = triangle_edge_rules(polynomial_order);
  std::fill(rhs, rhs + mesh.num_elements() * mode_count, 0.0);

  std::vector<double> basis(mode_count), gradient_xi(mode_count), gradient_eta(mode_count);
  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const TriangleMap map = element_map(mesh, element);
    const double *modes = coefficients + element * mode_count;
    double *element_rhs = rhs + element * mode_count;
    for (int q = 0; q < volume_rule.size(); ++q) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(volume_rule.xi[q], volume_rule.eta[q], x, y);
      const ScalarAdvection2D equation = velocity(x, y);
      dubiner_basis(volume_rule.xi[q], volume_rule.eta[q], polynomial_order, basis.data());
      dubiner_basis_gradient(volume_rule.xi[q], volume_rule.eta[q], polynomial_order,
                             gradient_xi.data(), gradient_eta.data());
      const double state = modal_sum(modes, basis.data(), mode_count);
      for (std::size_t mode = 0; mode < mode_count; ++mode) {
        double gradient_x = 0.0, gradient_y = 0.0;
        map.scaled_gradient(gradient_xi[mode], gradient_eta[mode], gradient_x, gradient_y);
        element_rhs[mode] += volume_rule.weights[q] * state *
                             (equation.velocity_x * gradient_x + equation.velocity_y * gradient_y);
      }
    }
  }

  std::vector<double> left_basis(mode_count), right_basis(mode_count);
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const TriangleFace &face = mesh.faces[face_index];
    const std::size_t left = static_cast<std::size_t>(face.left_element);
    const TriangleMap left_map = element_map(mesh, left);
    const TriangleEdgeRule &left_rule = edge_rules[static_cast<std::size_t>(face.left_local_edge)];
    const double edge_scale = face_edge_scale(left_map, face.left_local_edge);

    for (int q = 0; q < edge_point_count; ++q) {
      double x = 0.0, y = 0.0;
      left_map.reference_to_physical(left_rule.xi[q], left_rule.eta[q], x, y);
      const ScalarAdvection2D equation = velocity(x, y);
      dubiner_basis(left_rule.xi[q], left_rule.eta[q], polynomial_order, left_basis.data());
      const double left_state =
          modal_sum(coefficients + left * mode_count, left_basis.data(), mode_count);

      double right_state = 0.0;
      if (!face.is_boundary()) {
        double right_xi = 0.0, right_eta = 0.0;
        neighbor_edge_point(face, left_rule.xi[q], left_rule.eta[q], right_xi, right_eta);
        dubiner_basis(right_xi, right_eta, polynomial_order, right_basis.data());
        right_state =
            modal_sum(coefficients + static_cast<std::size_t>(face.right_element) * mode_count,
                      right_basis.data(), mode_count);
      } else {
        right_state = boundary_exterior_state[face_index * edge_point_count + q];
      }

      const double flux = upwind_normal_flux(equation, face.outward_normal_x, face.outward_normal_y,
                                             left_state, right_state);
      const double weight = left_rule.weights[q] * edge_scale;
      double *left_rhs = rhs + left * mode_count;
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        left_rhs[mode] -= weight * left_basis[mode] * flux;
      if (!face.is_boundary()) {
        double *right_rhs = rhs + static_cast<std::size_t>(face.right_element) * mode_count;
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          right_rhs[mode] += weight * right_basis[mode] * flux;
      }
    }
  }

  for (std::size_t element = 0; element < mesh.num_elements(); ++element) {
    const double inverse_mass = 1.0 / element_map(mesh, element).jacobian_determinant();
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      rhs[element * mode_count + mode] *= inverse_mass;
  }
}

inline void triangle_scalar_advection_residual(const TriangleMesh &mesh,
                                               const ScalarAdvection2D &equation,
                                               int polynomial_order, const double *coefficients,
                                               const double *boundary_exterior_state, double *rhs) {
  triangle_scalar_advection_residual(
      mesh, [equation](double, double) { return equation; }, polynomial_order, coefficients,
      boundary_exterior_state, rhs);
}

} // namespace wasl
