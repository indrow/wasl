// Physical tags of boundary edges in an MSH document.

#pragma once

#include <unordered_map>

#include "wasl/mesh/msh/msh_document.hpp"
#include "wasl/mesh/triangle/face_key.hpp"

namespace wasl {

class BoundaryEdgeTags {
public:
  explicit BoundaryEdgeTags(const MshDocument &document);

  // Tag of the boundary edge with this sorted key, or invalid_physical.
  int at(const FaceKey &key) const;

private:
  const MshDocument &document_;
  std::unordered_map<FaceKey, int, FaceKeyHash> line_element_tags_;
};

} // namespace wasl
