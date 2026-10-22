// Unpack the physical names stored in a prepared mesh.

#pragma once

#include "wasl/mesh/physical_groups.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// Requires a validated prepared mesh.
PhysicalNames prepared_physical_names(const PreparedMesh &mesh);

} // namespace wasl
