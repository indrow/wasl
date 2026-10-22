// A boundary condition attached to the boundary faces of one physical group.

#pragma once

#include <string>
#include <vector>

#include "wasl/boundaries/line_boundary_conditions.hpp"

namespace wasl {

struct NamedBoundary {
  std::string name;
  BoundaryCondition condition;
};

// An empty list means every boundary is periodic.
using NamedBoundaries = std::vector<NamedBoundary>;

} // namespace wasl
