// Zero-gradient outflow on both ends of a 1D scalar-advection mesh. See
// line_boundary_conditions.hpp for mixing policies per end.

#pragma once

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/equations/euler/euler.hpp"

namespace wasl {

inline EulerConserved euler_outflow_state(const EulerConserved &interior) noexcept {
  return interior;
}

inline ExteriorFunction make_outflow_exterior(const LineMesh &mesh, int polynomial_order) {
  const BoundaryCondition outflow{BoundaryKind::outflow, 0.0};
  return make_line_boundary_exterior(mesh, polynomial_order, {outflow, outflow});
}

} // namespace wasl
