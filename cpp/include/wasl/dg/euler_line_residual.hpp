// CPU modal DG residual for ideal-gas Euler on affine line elements.
// Coefficients use U[variable][mode][element]. Boundary traces use
// U[variable][face].

#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/reference_element.hpp"
#include "wasl/equations/euler/flux.hpp"
#include "wasl/mesh/line_geometry.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/quadrature/quadrature.hpp"

namespace wasl {
namespace detail {

inline double euler_line_component(const EulerConserved &state, int variable) noexcept {
  if (variable == 0)
    return state.density;
  if (variable == 1)
    return state.momentum_x;
  if (variable == 2)
    return state.momentum_y;
  return state.total_energy;
}

inline EulerConserved evaluate_euler_line_state(const double *coefficients,
                                                std::size_t element_count, std::size_t mode_count,
                                                std::size_t element, const double *basis) noexcept {
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
inline void euler_line_residual(const LineMesh &mesh, const Euler &equation,
                                NumericalFlux numerical_flux, int polynomial_order,
                                const double *coefficients, const double *boundary_exterior_state,
                                double *rhs) {
  const std::size_t element_count = mesh.num_elements();
  const std::size_t mode_count = line_mode_count(polynomial_order);
  const auto volume_rule = gauss_legendre(std::max(1, polynomial_order + 1));
  std::fill(rhs, rhs + Euler::variable_count * mode_count * element_count, 0.0);

  std::vector<double> basis(mode_count), gradient(mode_count);
  for (std::size_t element = 0; element < element_count; ++element) {
    for (int point = 0; point < volume_rule.size(); ++point) {
      legendre_basis(volume_rule.points[point], polynomial_order, basis.data());
      legendre_basis_gradient(volume_rule.points[point], polynomial_order, gradient.data());
      const auto state = detail::evaluate_euler_line_state(coefficients, element_count, mode_count,
                                                           element, basis.data());
      const auto flux = euler_flux_x(equation, state);
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + element] +=
              volume_rule.weights[point] * detail::euler_line_component(flux, variable) *
              gradient[mode];
    }
  }

  std::vector<double> left_basis(mode_count), right_basis(mode_count);
  legendre_basis(-1.0, polynomial_order, left_basis.data());
  legendre_basis(1.0, polynomial_order, right_basis.data());
  for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
    const auto &face = mesh.faces[face_index];
    const auto left = static_cast<std::size_t>(face.left_element);
    const double *owner_basis = face.outward_normal < 0.0 ? left_basis.data() : right_basis.data();
    const auto owner = detail::evaluate_euler_line_state(coefficients, element_count, mode_count,
                                                         left, owner_basis);
    EulerConserved exterior{};
    if (face.right_element != invalid_element) {
      exterior = detail::evaluate_euler_line_state(coefficients, element_count, mode_count,
                                                   static_cast<std::size_t>(face.right_element),
                                                   left_basis.data());
    } else {
      double *values[] = {&exterior.density, &exterior.momentum_x, &exterior.momentum_y,
                          &exterior.total_energy};
      for (int variable = 0; variable < Euler::variable_count; ++variable)
        *values[variable] =
            boundary_exterior_state[static_cast<std::size_t>(variable) * mesh.faces.size() +
                                    face_index];
    }
    const auto flux = face.outward_normal < 0.0
                          ? numerical_flux(equation, exterior, owner, 1.0, 0.0)
                          : numerical_flux(equation, owner, exterior, 1.0, 0.0);
    for (int variable = 0; variable < Euler::variable_count; ++variable) {
      const auto value = detail::euler_line_component(flux, variable);
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + left] +=
            (face.outward_normal < 0.0 ? 1.0 : -1.0) * owner_basis[mode] * value;
      if (face.right_element != invalid_element) {
        const auto right = static_cast<std::size_t>(face.right_element);
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + right] +=
              left_basis[mode] * value;
      }
    }
  }

  for (std::size_t element = 0; element < element_count; ++element) {
    const double inverse_jacobian = 1.0 / element_map(mesh, element).jacobian();
    for (int variable = 0; variable < Euler::variable_count; ++variable)
      for (std::size_t mode = 0; mode < mode_count; ++mode)
        rhs[(static_cast<std::size_t>(variable) * mode_count + mode) * element_count + element] *=
            inverse_jacobian;
  }
}

} // namespace wasl
