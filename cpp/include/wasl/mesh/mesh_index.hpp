// Sentinel values shared by every mesh topology.

#pragma once

#include <cstdint>

namespace wasl {

// The missing right element of a boundary face.
constexpr std::int32_t invalid_element = -1;

// An entity without a Gmsh physical group.
constexpr int invalid_physical = -1;

} // namespace wasl
