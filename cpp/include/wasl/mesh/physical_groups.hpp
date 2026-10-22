// Gmsh physical-group bookkeeping shared by line and triangle meshes.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace wasl {

// (dimension, tag) -> physical name.
using PhysicalNames = std::map<std::pair<int, int>, std::string>;

// Physical tag -> local element or face indices, ascending.
using PhysicalMembers = std::unordered_map<int, std::vector<std::int32_t>>;

} // namespace wasl
