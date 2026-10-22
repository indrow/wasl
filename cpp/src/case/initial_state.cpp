#include "wasl/case/initial_state.hpp"

#include <algorithm>

#include "case/initial_regions.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/case/scalar_coefficient_layout.hpp"
#include "wasl/dg/l2_projection.hpp"
#include "wasl/equations/euler/primitive.hpp"
#include "wasl/quadrature/gauss_legendre.hpp"
#include "wasl/quadrature/triangle_volume_rule.hpp"

namespace wasl {
namespace {

// Initial data may be discontinuous, so projections over-integrate: 7 Gauss
// points on lines and at least an order-5 rule on triangles.
QuadratureRule line_initial_rule() { return gauss_legendre(7); }

TriangleQuadratureRule triangle_initial_rule(int polynomial_order) {
  return triangle_volume_rule(std::max(polynomial_order, 5));
}

void store_conserved(const EulerConserved &state, double *values) {
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    values[variable] = euler_component(state, variable);
}

} // namespace

std::vector<double> initial_coefficients(const LineMesh &mesh, const ScalarAdvectionCase &config) {
  const auto projected = project_line_field(
      mesh, config.polynomial_order, line_initial_rule(), 1,
      [&](double x, double *values) { values[0] = scalar_initial_value(config, x); });
  return solver_coefficients(projected, mesh.num_elements(),
                             line_mode_count(config.polynomial_order));
}

std::vector<double> initial_coefficients(const TriangleMesh &mesh,
                                         const ScalarAdvectionCase &config) {
  const auto projected = project_triangle_field(
      mesh, config.polynomial_order, triangle_initial_rule(config.polynomial_order), 1,
      [&](double x, double y, double *values) { values[0] = scalar_initial_value(config, x, y); });
  return solver_coefficients(projected, mesh.num_elements(),
                             triangle_mode_count(config.polynomial_order));
}

std::vector<double> initial_coefficients(const LineMesh &mesh, const EulerCase &config,
                                         const Euler &equation) {
  return project_line_field(mesh, config.polynomial_order, line_initial_rule(),
                            Euler::variable_count, [&](double x, double *values) {
                              store_conserved(conserved(equation, euler_initial_value(config, x)),
                                              values);
                            });
}

std::vector<double> initial_coefficients(const TriangleMesh &mesh, const EulerCase &config,
                                         const Euler &equation) {
  return project_triangle_field(
      mesh, config.polynomial_order, triangle_initial_rule(config.polynomial_order),
      Euler::variable_count, [&](double x, double y, double *values) {
        store_conserved(conserved(equation, euler_initial_value(config, x, y)), values);
      });
}

} // namespace wasl
