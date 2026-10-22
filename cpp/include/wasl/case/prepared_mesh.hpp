#pragma once

#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {

// Rebuild the host mesh types used by residual kernels from a prepared mesh.
// Requires a single-rank prepared mesh (partition_rank_count == 1).
LineMesh line_mesh_from_prepared(const PreparedMesh &prepared);
TriangleMesh triangle_mesh_from_prepared(const PreparedMesh &prepared);

} // namespace wasl
