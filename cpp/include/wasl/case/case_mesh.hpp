// The prepared mesh of a case, with its canonical path for checkpoint
// provenance.

#pragma once

#include <filesystem>

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

struct CaseMesh {
  std::filesystem::path path; // weakly canonical
  PreparedMesh prepared;
};

// Throws when the mesh cannot be read or its dimension differs from the case.
CaseMesh load_case_mesh(const std::filesystem::path &prepared_mesh, int dimensions);

} // namespace wasl
