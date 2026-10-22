// Resolve named case boundaries to the two ends of a line mesh.

#pragma once

#include <vector>

#include "wasl/boundaries/line_boundary_conditions.hpp"
#include "wasl/case/named_boundary.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// Every boundary face's physical name must name a configured boundary, every
// configured boundary must match a boundary face, and periodic ends come in
// pairs. An empty configuration means periodic. Throws naming the offending
// boundary.
LineBoundaryConditions resolve_line_boundaries(const PreparedMesh &prepared, const LineMesh &mesh,
                                               const NamedBoundaries &boundaries);

// True when every configured boundary is periodic.
bool all_periodic(const NamedBoundaries &boundaries);

} // namespace wasl
