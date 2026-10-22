// Flatten a canonical mesh into a single-rank prepared mesh.

#pragma once

#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

PreparedMesh prepare_mesh(const LineMesh &mesh);
PreparedMesh prepare_mesh(const TriangleMesh &mesh);

} // namespace wasl
