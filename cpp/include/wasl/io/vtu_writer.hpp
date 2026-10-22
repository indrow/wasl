// ASCII VTK unstructured-grid output of cell values, without linking VTK.

#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// cell_values is packed [cell][component] and written as the "state" array.
void write_vtu(const std::filesystem::path &path, const PreparedMesh &mesh,
               const std::vector<double> &cell_values, std::size_t component_count = 1);

} // namespace wasl
