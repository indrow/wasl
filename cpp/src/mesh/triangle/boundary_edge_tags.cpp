#include "mesh/triangle/boundary_edge_tags.hpp"

#include <cstdint>
#include <cstdlib>

#include "wasl/mesh/mesh_index.hpp"

namespace wasl {

// Boundary line elements map each meshed edge to the physical tag of the curve
// that owns it. A curve split into many segments only lists its end nodes as
// bounding points, so line elements are authoritative whenever present.
BoundaryEdgeTags::BoundaryEdgeTags(const MshDocument &document) : document_(document) {
  for (const auto &element : document.elements) {
    if (element.element_type != msh_line_element || element.nodes.size() != 2)
      continue;
    line_element_tags_.emplace(make_face_key(element.nodes[0], element.nodes[1]),
                               document.first_physical_tag(element.entity_dim, element.entity_tag));
  }
}

// Without a line element, fall back to a single-segment curve whose two bounding
// points match the edge. Gmsh signs bounding points by orientation; only the
// magnitude is the node tag.
int BoundaryEdgeTags::at(const FaceKey &key) const {
  const auto found = line_element_tags_.find(key);
  if (found != line_element_tags_.end())
    return found->second;
  for (const auto &curve : document_.curve_entities) {
    if (curve.bounding_points.size() != 2)
      continue;
    const auto start = static_cast<std::int64_t>(std::abs(curve.bounding_points[0]));
    const auto end = static_cast<std::int64_t>(std::abs(curve.bounding_points[1]));
    if (make_face_key(start, end) == key)
      return curve.physical_tags.empty() ? invalid_physical : curve.physical_tags.front();
  }
  return invalid_physical;
}

} // namespace wasl
