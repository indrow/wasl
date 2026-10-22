// Structural validation of a prepared mesh at file and process boundaries.

#pragma once

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// Throws std::runtime_error naming the first inconsistent array or index.
void validate_prepared_mesh(const PreparedMesh &mesh);

} // namespace wasl
