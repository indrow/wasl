// Canonical triangle mesh topology for straight-sided 2D elements.
//
// Mirrors line_mesh.hpp: elements carry their physical vertices and node tags;
// faces are the canonical edges with deterministic, MSH-record-order-independent
// left/right ownership. Element connectivity is normalized to positive physical
// orientation (NUMERICAL_CONVENTIONS.md). Each face records, per side, the
// ordered pair of global vertex IDs as the adjacent element traverses its local
// edge, so a later residual can map neighbor traces into the canonical face
// coordinate.
//
// Local edges of a triangle element follow the reference convention
// (triangle_geometry.hpp):
//   edge 0: v0 -> v1,  edge 1: v1 -> v2,  edge 2: v2 -> v0.
// A face key is the sorted pair of its two global vertex IDs. Interior faces
// have exactly two adjacent elements; boundary faces have one owner and an
// invalid right element. The canonical normal points from the left element to
// the right element; boundary normals point outward.

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "wasl/mesh/line_mesh.hpp" // invalid_element, invalid_physical
#include "wasl/mesh/triangle_geometry.hpp"

namespace wasl {

struct MshDocument;

// A face (edge) key: the sorted pair of its two global vertex IDs. Equality and
// hashing are order-independent so an edge shared by two elements maps to one
// slot regardless of which element recorded it first.
struct FaceKey {
  std::int64_t first;  // smaller vertex ID
  std::int64_t second; // larger vertex ID

  bool operator==(const FaceKey &other) const noexcept {
    return first == other.first && second == other.second;
  }
};

struct FaceKeyHash {
  std::size_t operator()(const FaceKey &key) const noexcept {
    // Fold the two sorted vertex IDs into one hash, order-independent.
    const auto h1 = std::hash<std::int64_t>{}(key.first);
    const auto h2 = std::hash<std::int64_t>{}(key.second);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

struct TriangleElement {
  std::int64_t node0;
  std::int64_t node1;
  std::int64_t node2;
  double x0, y0; // image of reference v0 = (-1, -1)
  double x1, y1; // image of reference v1 = ( 1, -1)
  double x2, y2; // image of reference v2 = (-1,  1)
  int physical_tag;

  // The node tag at local vertex v, after orientation normalization.
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
};

struct TriangleFace {
  std::int64_t vertex_a; // sorted face key, first vertex ID
  std::int64_t vertex_b; // sorted face key, second vertex ID
  std::int32_t left_element;
  std::int32_t right_element; // invalid_element on a boundary face
  int left_local_edge;        // local edge (0,1,2) the left element uses
  int right_local_edge;       // local edge (0,1,2) the right element uses
  // Ordered vertex pair as each side traverses its local edge, for trace
  // orientation in the canonical face coordinate.
  std::int64_t left_vertex0;
  std::int64_t left_vertex1;
  std::int64_t right_vertex0;
  std::int64_t right_vertex1;
  double outward_normal_x; // canonical (interior) or outward (boundary) normal
  double outward_normal_y;
  int physical_tag; // boundary physical tag, or invalid_physical
};

struct TriangleMesh {
  std::vector<TriangleElement> elements;
  std::vector<TriangleFace> faces;
  std::map<std::pair<int, int>, std::string> physical_names;
  std::unordered_map<int, std::vector<std::int32_t>> physical_to_elements;
  std::unordered_map<int, std::vector<std::int32_t>> physical_to_faces;

  std::size_t num_elements() const noexcept { return elements.size(); }
};

TriangleMesh build_triangle_mesh(const MshDocument &document);
TriangleMap triangle_element_map(const TriangleMesh &mesh, std::size_t element);
std::size_t count_interior_faces(const TriangleMesh &mesh) noexcept;
std::size_t count_boundary_faces(const TriangleMesh &mesh) noexcept;

} // namespace wasl
