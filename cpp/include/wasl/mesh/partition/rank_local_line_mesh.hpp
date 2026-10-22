// A rank's local view of a global line mesh.

#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "wasl/mesh/line/line_mesh.hpp"

namespace wasl {

// Local elements are the owned elements followed by the ghost elements, so the
// residual reads a ghost neighbor exactly as an owned one. Only faces whose
// elements are all local are kept, renumbered to local indices.
struct RankLocalLineMesh {
  LineMesh mesh;
  std::map<std::int32_t, std::int32_t> global_to_local;
};

RankLocalLineMesh rank_local_line_mesh(const LineMesh &global,
                                       const std::vector<std::int32_t> &owned,
                                       const std::vector<std::int32_t> &ghosts);

} // namespace wasl
