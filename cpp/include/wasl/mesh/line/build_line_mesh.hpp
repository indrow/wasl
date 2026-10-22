// Build the canonical line mesh from an MSH document.

#pragma once

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/msh/msh_document.hpp"

namespace wasl {

// Orient every 2-node line element left to right, sort by position, and build
// endpoint faces. Throws for missing nodes, zero-length elements, or gaps.
LineMesh build_line_mesh(const MshDocument &document);

} // namespace wasl
