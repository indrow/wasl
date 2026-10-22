// The elements one rank owns, as a standalone prepared mesh for rank-local
// visualization pieces.

#pragma once

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// The slice has no faces and reports a single-rank partition. Throws when the
// rank is out of range or owns no elements.
PreparedMesh owned_mesh_slice(const PreparedMesh &mesh, int rank);

} // namespace wasl
