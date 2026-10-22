// Per-end boundary conditions for 1D scalar advection.
//
// A line mesh has exactly two boundary faces: the left end (outward normal -1)
// and the right end (outward normal +1). Each end is periodic, inflow, or
// outflow:
//   periodic  wraps to the opposite end's trace; both ends must be periodic.
//   inflow    prescribes a fixed exterior state. Where the wind blows outward
//             the upwind flux selects the interior trace, so the value is
//             unused there.
//   outflow   copies the owner's own trace (zero gradient), so waves leave
//             freely.
// The exterior is recomputed from each stage's coefficients, so boundary values
// are applied at every SSP-RK3 stage.

#pragma once

#include <cstddef>
#include <vector>

#include "wasl/boundaries/line_exterior.hpp"

namespace wasl {

enum class BoundaryKind { periodic, inflow, outflow };

struct BoundaryCondition {
  BoundaryKind kind = BoundaryKind::periodic;
  double value = 0.0; // inflow state; unused otherwise
};

struct LineBoundaryConditions {
  BoundaryCondition left;
  BoundaryCondition right;
};

inline ExteriorFunction make_line_boundary_exterior(const LineMesh &mesh, int polynomial_order,
                                                    LineBoundaryConditions conditions) {
  return make_line_exterior(
      mesh, polynomial_order,
      [&mesh, polynomial_order, conditions](const std::vector<double> &coefficients,
                                            const LineFace &face) {
        const bool left_end = face.outward_normal < 0.0;
        const BoundaryCondition &condition = left_end ? conditions.left : conditions.right;
        switch (condition.kind) {
        case BoundaryKind::inflow:
          return condition.value;
        case BoundaryKind::outflow:
          return line_element_trace(coefficients, static_cast<std::size_t>(face.left_element),
                                    polynomial_order, left_end ? -1.0 : 1.0);
        case BoundaryKind::periodic:
        default:
          return left_end ? line_element_trace(coefficients, mesh.num_elements() - 1,
                                               polynomial_order, 1.0)
                          : line_element_trace(coefficients, 0, polynomial_order, -1.0);
        }
      });
}

} // namespace wasl
