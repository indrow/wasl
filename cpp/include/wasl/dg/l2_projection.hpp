// L2 projection of a physical field onto the orthonormal modal basis.
//
// With an orthonormal basis the reference mass matrix is the identity, so each
// coefficient is the quadrature sum c_i = sum_q w_q f(x(q)) phi_i(q). Results
// use U[variable][mode][element].

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/basis/dubiner_triangle.hpp"
#include "wasl/basis/legendre.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/dg/modal_layout.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"
#include "wasl/quadrature/quadrature_rule.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {
namespace detail {

inline void accumulate_projection(const ModalLayout &layout, std::size_t variable_count,
                                  std::size_t element, double weight, const double *values,
                                  const double *basis, double *coefficients) noexcept {
  for (std::size_t variable = 0; variable < variable_count; ++variable)
    for (std::size_t mode = 0; mode < layout.mode_count; ++mode)
      coefficients[layout.index(variable, mode, element)] +=
          weight * values[variable] * basis[mode];
}

} // namespace detail

// Field: (x, double *values) filling variable_count values.
template <typename Field>
std::vector<double> project_line_field(const LineMesh &mesh, int polynomial_order,
                                       const QuadratureRule &quadrature, std::size_t variable_count,
                                       Field field) {
  const ModalLayout layout{line_mode_count(polynomial_order), mesh.num_elements()};
  std::vector<double> coefficients(variable_count * layout.mode_count * layout.element_count, 0.0);
  std::vector<double> basis(layout.mode_count), values(variable_count);
  for (std::size_t element = 0; element < layout.element_count; ++element) {
    const auto map = element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      const auto xi = quadrature.points[point];
      legendre_basis(xi, polynomial_order, basis.data());
      field(map.reference_to_physical(xi), values.data());
      detail::accumulate_projection(layout, variable_count, element, quadrature.weights[point],
                                    values.data(), basis.data(), coefficients.data());
    }
  }
  return coefficients;
}

// Field: (x, y, double *values) filling variable_count values.
template <typename Field>
std::vector<double> project_triangle_field(const TriangleMesh &mesh, int polynomial_order,
                                           const TriangleQuadratureRule &quadrature,
                                           std::size_t variable_count, Field field) {
  const ModalLayout layout{triangle_mode_count(polynomial_order), mesh.num_elements()};
  std::vector<double> coefficients(variable_count * layout.mode_count * layout.element_count, 0.0);
  std::vector<double> basis(layout.mode_count), values(variable_count);
  for (std::size_t element = 0; element < layout.element_count; ++element) {
    const auto map = element_map(mesh, element);
    for (int point = 0; point < quadrature.size(); ++point) {
      double x = 0.0, y = 0.0;
      map.reference_to_physical(quadrature.xi[point], quadrature.eta[point], x, y);
      dubiner_basis(quadrature.xi[point], quadrature.eta[point], polynomial_order, basis.data());
      field(x, y, values.data());
      detail::accumulate_projection(layout, variable_count, element, quadrature.weights[point],
                                    values.data(), basis.data(), coefficients.data());
    }
  }
  return coefficients;
}

} // namespace wasl
