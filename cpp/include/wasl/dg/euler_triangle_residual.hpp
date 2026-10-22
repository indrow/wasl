// CPU modal DG residual for two-dimensional Euler on affine triangles.
// Coefficients use U[variable][mode][element]. Boundary traces use
// U[variable][face][edge quadrature point].

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/flux.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/quadrature/triangle.hpp"

namespace wasl {
namespace detail {

inline double euler_component(const EulerConserved &state, int variable) noexcept {
  if (variable == 0)
    return state.density;
  if (variable == 1)
    return state.momentum_x;
  if (variable == 2)
    return state.momentum_y;
  return state.total_energy;
}

inline EulerConserved evaluate_euler_state(const double *coefficients, std::size_t element_count,
                                           std::size_t mode_count, std::size_t element,
                                           const double *basis) noexcept {
  EulerConserved state{0.0, 0.0, 0.0, 0.0};
  double *values[] = {&state.density, &state.momentum_x, &state.momentum_y, &state.total_energy};
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      *values[variable] +=
          coefficients[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count +
                       element] *
          basis[mode];
  return state;
}

} // namespace detail

template <typename NumericalFlux>
inline void euler_triangle_residual(const TriangleMesh &mesh, const Euler &equation,
                                    NumericalFlux numerical_flux, int polynomial_order,
                                    const double *coefficients,
                                    const double *boundary_exterior_state, double *rhs) {
  const std::size_t element_count = mesh.num_elements();
  const std::size_t mode_count = triangle_mode_count(polynomial_order);
  const int edge_point_count = polynomial_order + 1;
  const auto volume_rule = triangle_volume_rule(polynomial_order + 2);
  const auto edge_rules = triangle_edge_rules(polynomial_order);
  std::fill(rhs, rhs + Euler::variable_count * mode_count * element_count, 0.0);

  std::vector<double> basis(mode_count), gradient_xi(mode_count), gradient_eta(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    const TriangleMap map = triangle_element_map(mesh, element);
    for (int q = 0; q < volume_rule.size(); ++q) {
      dubiner_basis(volume_rule.xi[q], volume_rule.eta[q], polynomial_order, basis.data());
      dubiner_basis_gradient(volume_rule.xi[q], volume_rule.eta[q], polynomial_order,
                             gradient_xi.data(), gradient_eta.data());
      const EulerConserved state = detail::evaluate_euler_state(coefficients, element_count,
                                                                mode_count, element, basis.data());
      const EulerConserved flux_x = euler_flux_x(equation, state);
      const EulerConserved flux_y = euler_flux_y(equation, state);
      for (std::size_t mode = 0; mode < mode_count; ++mode) {
        const double scaled_gradient_x =
            0.5 * (map.by() * gradient_xi[mode] - map.ay() * gradient_eta[mode]);
        const double scaled_gradient_y =
            0.5 * (-map.bx() * gradient_xi[mode] + map.ax() * gradient_eta[mode]);
        for (int variable = 0; variable < Euler::variable_count; ++variable) {
          rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + element] +=
              volume_rule.weights[q] *
              (detail::euler_component(flux_x, variable) * scaled_gradient_x +
               detail::euler_component(flux_y, variable) * scaled_gradient_y);
        }
      }
    }
  }

  std::vector<double> left_basis(mode_count), right_basis(mode_count);
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const TriangleFace &face = mesh.faces[face_index];
    const std::size_t left = static_cast<std::size_t>(face.left_element);
    const TriangleMap left_map = triangle_element_map(mesh, left);
    const TriangleEdgeRule &left_rule = edge_rules[static_cast<std::size_t>(face.left_local_edge)];
    const double reference_length = face.left_local_edge == 1 ? 2.0 * std::sqrt(2.0) : 2.0;
    const double edge_scale = left_map.edge_length(face.left_local_edge) / reference_length;

    for (int q = 0; q < edge_point_count; ++q) {
      dubiner_basis(left_rule.xi[q], left_rule.eta[q], polynomial_order, left_basis.data());
      const EulerConserved left_state = detail::evaluate_euler_state(
          coefficients, element_count, mode_count, left, left_basis.data());
      EulerConserved right_state{};
      if (face.right_element != invalid_element) {
        const std::size_t right = static_cast<std::size_t>(face.right_element);
        double parameter = face.left_local_edge == 0   ? left_rule.xi[q]
                           : face.left_local_edge == 1 ? left_rule.eta[q]
                                                       : -left_rule.eta[q];
        if (face.left_vertex0 != face.right_vertex0)
          parameter = -parameter;
        double right_xi = -1.0;
        double right_eta = -1.0;
        if (face.right_local_edge == 0) {
          right_xi = parameter;
        } else if (face.right_local_edge == 1) {
          right_xi = -parameter;
          right_eta = parameter;
        } else {
          right_eta = -parameter;
        }
        dubiner_basis(right_xi, right_eta, polynomial_order, right_basis.data());
        right_state = detail::evaluate_euler_state(coefficients, element_count, mode_count, right,
                                                   right_basis.data());
      } else {
        double *values[] = {&right_state.density, &right_state.momentum_x, &right_state.momentum_y,
                            &right_state.total_energy};
        for (int variable = 0; variable < Euler::variable_count; ++variable)
          *values[variable] =
              boundary_exterior_state[(static_cast<std::size_t>(variable) * mesh.faces.size() +
                                       face_index) *
                                          edge_point_count +
                                      q];
      }

      const EulerConserved flux = numerical_flux(equation, left_state, right_state,
                                                 face.outward_normal_x, face.outward_normal_y);
      const double weight = left_rule.weights[q] * edge_scale;
      for (int variable = 0; variable < Euler::variable_count; ++variable) {
        const double value = detail::euler_component(flux, variable) * weight;
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + left] -=
              value * left_basis[mode];
        if (face.right_element != invalid_element) {
          const std::size_t right = static_cast<std::size_t>(face.right_element);
          for (std::size_t mode = 0; mode < mode_count; ++mode)
            rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + right] +=
                value * right_basis[mode];
        }
      }
    }
  }

  for (std::size_t element = 0; element < element_count; ++element) {
    const double inverse_jacobian =
        1.0 / triangle_element_map(mesh, element).jacobian_determinant();
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + element] *=
            inverse_jacobian;
  }
}

} // namespace wasl
