// Exterior Euler states at the two ends of a line mesh.

#pragma once

#include <vector>

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/mesh/line/line_mesh.hpp"

namespace wasl {

// Exterior per [variable][face] from U[variable][mode][element]; interior
// entries stay zero. A periodic left end wraps to the last element's trace at
// xi = +1 and a periodic right end to the first element's trace at xi = -1. An
// outflow end copies its owner's trace, so every characteristic leaves through
// the upwind flux without reflection. Euler inflow throws: a prescribed state
// needs all four conserved variables.
std::vector<double> line_euler_exterior(const LineMesh &mesh, int polynomial_order,
                                        const LineBoundaryConditions &conditions,
                                        const std::vector<double> &coefficients);

} // namespace wasl
