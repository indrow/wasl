// Modal DG residual for ideal-gas Euler on affine line elements.
//
// Coefficients and rhs use U[variable][mode][element]; boundary exterior states
// use U[variable][face]. The weak form matches the scalar line residual: the
// returned rhs already includes the orthonormal mass inverse 1/J.

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/equations/euler/flux.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/quadrature/gauss_legendre.hpp"

namespace wasl {

// NumericalFlux: (Euler, left, right, normal_x, normal_y) -> EulerConserved.
template <typename NumericalFlux>
inline void euler_line_residual(const LineMesh &mesh, const Euler &equation,
                                NumericalFlux numerical_flux, int polynomial_order,
                                const double *coefficients, const double *boundary_exterior_state,
                                double *rhs) {
  const ModalLayout layout{line_mode_count(polynomial_order), mesh.num_elements()};
  const auto volume_rule = gauss_legendre(std::max(1, polynomial_order + 1));
  std::fill(rhs, rhs + Euler::variable_count * layout.mode_count * layout.element_count, 0.0);

  std::vector<double> basis(layout.mode_count), gradient(layout.mode_count);
  for (std::size_t element = 0; element < layout.element_count; ++element) {
    for (int point = 0; point < volume_rule.size(); ++point) {
      legendre_basis(volume_rule.points[point], polynomial_order, basis.data());
      legendre_basis_gradient(volume_rule.points[point], polynomial_order, gradient.data());
      const auto flux =
          euler_flux_x(equation, evaluate_euler_state(layout, coefficients, element, basis.data()));
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
          rhs[layout.index(static_cast<std::size_t>(variable), mode, element)] +=
              volume_rule.weights[point] * euler_component(flux, variable) * gradient[mode];
    }
  }

  std::vector<double> left_basis(layout.mode_count), right_basis(layout.mode_count);
  legendre_basis(-1.0, polynomial_order, left_basis.data());
  legendre_basis(1.0, polynomial_order, right_basis.data());
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    const auto left = static_cast<std::size_t>(face.left_element);
    // Only the left boundary face (normal -1) is owned at the owner's xi = -1.
    const bool owned_at_left_end = face.outward_normal < 0.0;
    const double *owner_basis = owned_at_left_end ? left_basis.data() : right_basis.data();
    const auto owner = evaluate_euler_state(layout, coefficients, left, owner_basis);
    EulerConserved exterior{};
    if (!face.is_boundary()) {
      exterior = evaluate_euler_state(
          layout, coefficients, static_cast<std::size_t>(face.right_element), left_basis.data());
    } else {
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        euler_component(exterior, variable) =
            boundary_exterior_state[static_cast<std::size_t>(variable) * mesh.faces.size() +
                                    face_index];
    }
    // The flux is always evaluated left-to-right along +x.
    const auto flux = owned_at_left_end ? numerical_flux(equation, exterior, owner, 1.0, 0.0)
                                        : numerical_flux(equation, owner, exterior, 1.0, 0.0);
    for (int variable = 0; variable < Euler::variable_count; ++variable) {
      const auto value = euler_component(flux, variable);
      for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
        rhs[layout.index(static_cast<std::size_t>(variable), mode, left)] +=
            (owned_at_left_end ? 1.0 : -1.0) * owner_basis[mode] * value;
      if (!face.is_boundary()) {
        const auto right = static_cast<std::size_t>(face.right_element);
        for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
          rhs[layout.index(static_cast<std::size_t>(variable), mode, right)] +=
              left_basis[mode] * value;
      }
    }
  }

  for (std::size_t element = 0; element < layout.element_count; ++element)
    scale_element(layout, Euler::variable_count, rhs, element,
                  1.0 / element_map(mesh, element).jacobian());
}

} // namespace wasl
