// Tag-indexed access to MSH nodes for the mesh builders.

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>

#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/physical_groups.hpp"

namespace wasl {

class MshNodeLookup {
public:
  // context prefixes error messages, for example "line mesh".
  MshNodeLookup(const MshDocument &document, std::string context);

  const MshNode &at(std::int64_t tag) const;

private:
  std::unordered_map<std::int64_t, const MshNode *> nodes_;
  std::string context_;
};

PhysicalNames physical_names_of(const MshDocument &document);

} // namespace wasl
