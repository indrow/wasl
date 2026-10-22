// Gmsh MSH 4.1 ASCII reader.

#pragma once

#include <string>

#include "wasl/mesh/msh/msh_document.hpp"

namespace wasl {

// Parse the sections Wasl uses and skip the rest. Throws std::runtime_error
// naming the offending field for malformed or unsupported input.
MshDocument read_msh(const std::string &path);

} // namespace wasl
