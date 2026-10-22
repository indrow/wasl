// Axis-aligned bounding box of a triangle mesh.

#pragma once

#include <algorithm>
#include <limits>

#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

struct BoundingBox {
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();

  double width() const noexcept { return max_x - min_x; }
  double height() const noexcept { return max_y - min_y; }
};

inline BoundingBox bounding_box(const TriangleMesh &mesh) noexcept {
  BoundingBox box;
  for (const auto &element : mesh.elements) {
    box.min_x = std::min({box.min_x, element.x0, element.x1, element.x2});
    box.max_x = std::max({box.max_x, element.x0, element.x1, element.x2});
    box.min_y = std::min({box.min_y, element.y0, element.y1, element.y2});
    box.max_y = std::max({box.max_y, element.y0, element.y1, element.y2});
  }
  return box;
}

} // namespace wasl
