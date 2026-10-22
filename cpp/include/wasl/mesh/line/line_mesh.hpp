// Canonical 1D line-mesh topology.
//
// Elements are sorted by x_left. Faces are the element endpoints: the left
// boundary face (normal -1), one interior face per adjacent pair (normal +1,
// left element first), and the right boundary face (normal +1). A boundary face
// stores invalid_element as its right element.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/mesh/line/line_geometry.hpp"
#include "wasl/mesh/mesh_index.hpp"
#include "wasl/mesh/physical_groups.hpp"

namespace wasl {

struct LineElement {
  std::int64_t left_node;
  std::int64_t right_node;
  double x_left;
  double x_right;
  int physical_tag;

  LineMap map() const noexcept { return LineMap{x_left, x_right}; }
};

struct LineFace {
  std::int32_t left_element;
  std::int32_t right_element;
  double x;
  double outward_normal;
  int physical_tag;

  bool is_boundary() const noexcept { return right_element == invalid_element; }
};

struct LineMesh {
  std::vector<LineElement> elements;
  std::vector<LineFace> faces;
  PhysicalNames physical_names;
  PhysicalMembers physical_to_elements;
  PhysicalMembers physical_to_faces;

  std::size_t num_elements() const noexcept { return elements.size(); }
};

LineMap element_map(const LineMesh &mesh, std::size_t element);
std::size_t count_interior_faces(const LineMesh &mesh) noexcept;
std::size_t count_boundary_faces(const LineMesh &mesh) noexcept;

} // namespace wasl
