#include "mesh/prepared/prepared_node_table.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>

namespace wasl {

std::int32_t PreparedNodeTable::add(std::int64_t tag, double x, double y) {
  const auto found = indices_.find(tag);
  if (found != indices_.end())
    return found->second;
  if (mesh_.node_tags.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    throw std::runtime_error("prepared mesh: local node count exceeds int32 range");
  const auto index = static_cast<std::int32_t>(mesh_.node_tags.size());
  indices_.emplace(tag, index);
  mesh_.node_tags.push_back(tag);
  mesh_.points.insert(mesh_.points.end(), {x, y, 0.0});
  return index;
}

} // namespace wasl
