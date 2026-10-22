// Modal DG residual for two-dimensional Euler on affine triangles.
//
// Coefficients and rhs use U[variable][mode][element]; boundary exterior states
// use U[variable][face][edge quadrature point]. The returned rhs includes the
// orthonormal mass inverse 1/det(J).

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/dg/triangle_face_trace.hpp"
#include "wasl/equations/euler/flux.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"
#include "wasl/quadrature/triangle_edge_rules.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {

// NumericalFlux: (Euler, left, right, normal_x, normal_y) -> EulerConserved.
template <typename NumericalFlux>
inline void euler_triangle_residual(const TriangleMesh &mesh, const Euler &equation,
                                    NumericalFlux numerical_flux, int polynomial_order,
                                    const double *coefficients,
                                    const double *boundary_exterior_state, double *rhs) {
  const ModalLayout layout{triangle_mode_count(polynomial_order), mesh.num_elements()};
  const int edge_point_count = polynomial_order + 1;
  const auto volume_rule = triangle_volume_rule(polynomial_order + 2);
  const auto edge_rules = triangle_edge_rules(polynomial_order);
  std::fill(rhs, rhs + Euler::variable_count * layout.mode_count * layout.element_count, 0.0);

  std::vector<double> basis(layout.mode_count), gradient_xi(layout.mode_count),
      gradient_eta(layout.mode_count);
  for (std::size_t element = 0; element < layout.element_count; ++element) {
    const TriangleMap map = element_map(mesh, element);
    for (int q = 0; q < volume_rule.size(); ++q) {
      dubiner_basis(volume_rule.xi[q], volume_rule.eta[q], polynomial_order, basis.data());
      dubiner_basis_gradient(volume_rule.xi[q], volume_rule.eta[q], polynomial_order,
                             gradient_xi.data(), gradient_eta.data());
      const EulerConserved state =
          evaluate_euler_state(layout, coefficients, element, basis.data());
      const EulerConserved flux_x = euler_flux_x(equation, state);
      const EulerConserved flux_y = euler_flux_y(equation, state);
      for (std::size_t mode = 0; mode < layout.mode_count; ++mode) {
        double gradient_x = 0.0, gradient_y = 0.0;
        map.scaled_gradient(gradient_xi[mode], gradient_eta[mode], gradient_x, gradient_y);
        for (int variable = 0; variable < Euler::variable_count; ++variable)
          rhs[layout.index(static_cast<std::size_t>(variable), mode, element)] +=
              volume_rule.weights[q] * (euler_component(flux_x, variable) * gradient_x +
                                        euler_component(flux_y, variable) * gradient_y);
      }
    }
  }

  std::vector<double> left_basis(layout.mode_count), right_basis(layout.mode_count);
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const TriangleFace &face = mesh.faces[face_index];
    const std::size_t left = static_cast<std::size_t>(face.left_element);
    const TriangleEdgeRule &left_rule = edge_rules[static_cast<std::size_t>(face.left_local_edge)];
    const double edge_scale = face_edge_scale(element_map(mesh, left), face.left_local_edge);

    for (int q = 0; q < edge_point_count; ++q) {
      dubiner_basis(left_rule.xi[q], left_rule.eta[q], polynomial_order, left_basis.data());
      const EulerConserved left_state =
          evaluate_euler_state(layout, coefficients, left, left_basis.data());
      EulerConserved right_state{};
      if (!face.is_boundary()) {
        double right_xi = 0.0, right_eta = 0.0;
        neighbor_edge_point(face, left_rule.xi[q], left_rule.eta[q], right_xi, right_eta);
        dubiner_basis(right_xi, right_eta, polynomial_order, right_basis.data());
        right_state = evaluate_euler_state(
            layout, coefficients, static_cast<std::size_t>(face.right_element), right_basis.data());
      } else {
        for (int variable = 0; variable < Euler::variable_count; ++variable)
          euler_component(right_state, variable) =
              boundary_exterior_state[(static_cast<std::size_t>(variable) * mesh.faces.size() +
                                       face_index) *
                                          edge_point_count +
                                      q];
      }

      const EulerConserved flux = numerical_flux(equation, left_state, right_state,
                                                 face.outward_normal_x, face.outward_normal_y);
      const double weight = left_rule.weights[q] * edge_scale;
      for (int variable = 0; variable < Euler::variable_count; ++variable) {
        const double value = euler_component(flux, variable) * weight;
        for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
          rhs[layout.index(static_cast<std::size_t>(variable), mode, left)] -=
              value * left_basis[mode];
        if (!face.is_boundary()) {
          const std::size_t right = static_cast<std::size_t>(face.right_element);
          for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
            rhs[layout.index(static_cast<std::size_t>(variable), mode, right)] +=
                value * right_basis[mode];
        }
      }
    }
  }

  for (std::size_t element = 0; element < layout.element_count; ++element)
    scale_element(layout, Euler::variable_count, rhs, element,
                  1.0 / element_map(mesh, element).jacobian_determinant());
}

} // namespace wasl
