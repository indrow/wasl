// Build the canonical triangle mesh from an MSH document.

#pragma once

#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

// Orient every 3-node triangle counterclockwise, build shared edges with
// deterministic ownership, and tag boundary edges from MSH line elements or
// curve entities. Throws for missing nodes, zero-area elements, or edges with
// more than two adjacent elements.
TriangleMesh build_triangle_mesh(const MshDocument &document);

} // namespace wasl
