// Shared assembly of 1D scalar exterior traces. Interior faces always use the
// right element's trace at xi = -1; only boundary faces differ by policy.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/legendre.hpp"
#include "wasl/basis/modal_expansion.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/exterior_function.hpp"
#include "wasl/mesh/line/line_mesh.hpp"

namespace wasl {

// Trace of one element (U[element][mode]) at reference endpoint xi.
inline double line_element_trace(const std::vector<double> &coefficients, std::size_t element,
                                 int polynomial_order, double xi) {
  const std::size_t mode_count = line_mode_count(polynomial_order);
  std::vector<double> basis(mode_count);
  legendre_basis(xi, polynomial_order, basis.data());
  return modal_sum(coefficients.data() + element * mode_count, basis.data(), mode_count);
}

// BoundaryValue: (coefficients, face) -> exterior value of a boundary face.
template <typename BoundaryValue>
ExteriorFunction make_line_exterior(const LineMesh &mesh, int polynomial_order,
                                    BoundaryValue boundary_value) {
  return [&mesh, polynomial_order,
          boundary_value](const std::vector<double> &coefficients) -> std::vector<double> {
    std::vector<double> exterior(mesh.faces.size(), 0.0);
    for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
      const auto &value = mesh.faces[face];
      exterior[face] =
          value.is_boundary()
              ? boundary_value(coefficients, value)
              : line_element_trace(coefficients, static_cast<std::size_t>(value.right_element),
                                   polynomial_order, -1.0);
    }
    return exterior;
  };
}

} // namespace wasl
