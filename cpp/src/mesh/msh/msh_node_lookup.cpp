#include "mesh/msh/msh_node_lookup.hpp"

#include <stdexcept>
#include <utility>

namespace wasl {

MshNodeLookup::MshNodeLookup(const MshDocument &document, std::string context)
    : context_(std::move(context)) {
  nodes_.reserve(document.nodes.size());
  for (const auto &node : document.nodes) {
    if (!nodes_.emplace(node.tag, &node).second)
      throw std::runtime_error(context_ + ": duplicate node tag " + std::to_string(node.tag));
  }
}

const MshNode &MshNodeLookup::at(std::int64_t tag) const {
  const auto found = nodes_.find(tag);
  if (found == nodes_.end())
    throw std::runtime_error(context_ + ": element references missing node " + std::to_string(tag));
  return *found->second;
}

PhysicalNames physical_names_of(const MshDocument &document) {
  PhysicalNames names;
  for (const auto &physical : document.physical_names)
    names[{physical.dim, physical.tag}] = physical.name;
  return names;
}

} // namespace wasl
