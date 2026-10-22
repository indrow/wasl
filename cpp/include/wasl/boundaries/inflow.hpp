// Inflow on both ends of a 1D scalar-advection mesh. See
// line_boundary_conditions.hpp for mixing policies per end.

#pragma once

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/equations/euler/euler.hpp"

namespace wasl {

inline EulerConserved euler_inflow_state(const EulerConserved &prescribed) noexcept {
  return prescribed;
}

inline ExteriorFunction make_inflow_exterior(const LineMesh &mesh, int polynomial_order,
                                             double inflow_value) {
  const BoundaryCondition inflow{BoundaryKind::inflow, inflow_value};
  return make_line_boundary_exterior(mesh, polynomial_order, {inflow, inflow});
}

} // namespace wasl
