// Versioned HDF5 prepared-mesh files.
//
// Format history:
//   1: connectivity, metrics, faces, physical names.
//   2: adds partition ownership (rank count, element and face owners).
// Version 1 reads back-fill a single-rank partition; writes always emit 2.

#pragma once

#include <filesystem>

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

void write_prepared_mesh(const std::filesystem::path &path, const PreparedMesh &mesh);
PreparedMesh read_prepared_mesh(const std::filesystem::path &path);

} // namespace wasl
