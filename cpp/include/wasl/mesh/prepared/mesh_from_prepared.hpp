// Rebuild the canonical meshes residual kernels use from a prepared mesh.

#pragma once

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// Both require a single-rank prepared mesh of the matching dimension.
LineMesh line_mesh_from_prepared(const PreparedMesh &prepared);
TriangleMesh triangle_mesh_from_prepared(const PreparedMesh &prepared);

} // namespace wasl
