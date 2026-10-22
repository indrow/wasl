// Deduplicating node insertion while assembling a prepared mesh.

#pragma once

#include <cstdint>
#include <map>

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

class PreparedNodeTable {
public:
  explicit PreparedNodeTable(PreparedMesh &mesh) : mesh_(mesh) {}

  // Local index of the node with this tag, appending it on first use.
  std::int32_t add(std::int64_t tag, double x, double y);

  // Local index of a node already added; throws std::out_of_range otherwise.
  std::int32_t at(std::int64_t tag) const { return indices_.at(tag); }

private:
  PreparedMesh &mesh_;
  std::map<std::int64_t, std::int32_t> indices_;
};

} // namespace wasl
