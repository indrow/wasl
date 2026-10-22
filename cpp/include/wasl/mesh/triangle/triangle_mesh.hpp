// Canonical triangle-mesh topology for straight-sided 2D elements.
//
// Elements carry their physical vertices and node tags, normalized to positive
// orientation. Local edges follow the reference convention
//   edge 0: v0 -> v1,  edge 1: v1 -> v2,  edge 2: v2 -> v0.
// Faces have deterministic, MSH-record-order-independent ownership: an interior
// face's left element is the smaller element index, and a boundary face stores
// invalid_element as its right element. The stored normal points from left to
// right, or outward on a boundary. Each face records the ordered vertex pair
// each side traverses so residuals can map neighbor traces into the canonical
// face coordinate.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "wasl/mesh/mesh_index.hpp"
#include "wasl/mesh/physical_groups.hpp"
#include "wasl/mesh/triangle/triangle_geometry.hpp"

namespace wasl {

// Local (start, end) vertices of each triangle edge.
constexpr int triangle_edge_vertices[3][2] = {{0, 1}, {1, 2}, {2, 0}};

struct TriangleElement {
  std::int64_t node0;
  std::int64_t node1;
  std::int64_t node2;
  double x0, y0; // image of reference v0 = (-1, -1)
  double x1, y1; // image of reference v1 = ( 1, -1)
  double x2, y2; // image of reference v2 = (-1,  1)
  int physical_tag;

  std::int64_t node(int local_vertex) const noexcept {
    switch (local_vertex) {
    case 0:
      return node0;
    case 1:
      return node1;
    default:
      return node2;
    }
  }

  TriangleMap map() const noexcept { return TriangleMap{x0, y0, x1, y1, x2, y2}; }
};

struct TriangleFace {
  std::int64_t vertex_a; // sorted face key, first vertex ID
  std::int64_t vertex_b; // sorted face key, second vertex ID
  std::int32_t left_element;
  std::int32_t right_element; // invalid_element on a boundary face
  int left_local_edge;
  int right_local_edge;
  // Ordered vertex pair as each side traverses its local edge.
  std::int64_t left_vertex0;
  std::int64_t left_vertex1;
  std::int64_t right_vertex0;
  std::int64_t right_vertex1;
  double outward_normal_x;
  double outward_normal_y;
  int physical_tag; // boundary physical tag, or invalid_physical

  bool is_boundary() const noexcept { return right_element == invalid_element; }
};

struct TriangleMesh {
  std::vector<TriangleElement> elements;
  std::vector<TriangleFace> faces;
  PhysicalNames physical_names;
  PhysicalMembers physical_to_elements;
  PhysicalMembers physical_to_faces;

  std::size_t num_elements() const noexcept { return elements.size(); }
};

TriangleMap element_map(const TriangleMesh &mesh, std::size_t element);
std::size_t count_interior_faces(const TriangleMesh &mesh) noexcept;
std::size_t count_boundary_faces(const TriangleMesh &mesh) noexcept;

} // namespace wasl
