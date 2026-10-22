// Periodic boundary policy for 1D scalar advection.
//
// The domain is a ring: the left boundary face wraps to the last element's
// trace at xi = +1, and the right boundary face wraps to the first element's
// trace at xi = -1. Total mass is conserved and a smooth pulse returns to its
// initial state after one revolution.

#pragma once

#include "wasl/boundaries/line_boundary_conditions.hpp"

namespace wasl {

inline ExteriorFunction make_periodic_exterior(const LineMesh &mesh, int polynomial_order) {
  return make_line_boundary_exterior(mesh, polynomial_order, {});
}

} // namespace wasl
