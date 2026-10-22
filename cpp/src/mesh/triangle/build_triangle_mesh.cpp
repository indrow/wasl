#include "wasl/mesh/triangle/build_triangle_mesh.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "mesh/msh/msh_node_lookup.hpp"
#include "mesh/triangle/boundary_edge_tags.hpp"
#include "wasl/mesh/triangle/face_key.hpp"

namespace wasl {
namespace {

// One element's view of an edge: which local edge it is and the vertex order
// the element traverses it in.
struct EdgeSide {
  std::int32_t element = invalid_element;
  int local_edge = 0;
  std::int64_t vertex0 = 0;
  std::int64_t vertex1 = 0;
};

struct EdgeAdjacency {
  FaceKey key{};
  EdgeSide first;
  EdgeSide second; // element stays invalid until a second element shares the edge
};

using EdgeTable = std::unordered_map<FaceKey, EdgeAdjacency, FaceKeyHash>;

// Triangles with counterclockwise vertices: a negative signed area swaps the
// trailing two vertices.
std::vector<TriangleElement> oriented_elements(const MshDocument &document,
                                               const MshNodeLookup &nodes) {
  std::vector<TriangleElement> elements;
  elements.reserve(document.elements.size());
  for (const auto &element : document.elements) {
    if (element.element_type != msh_triangle_element || element.nodes.size() != 3)
      continue;
    std::array<const MshNode *, 3> v = {&nodes.at(element.nodes[0]), &nodes.at(element.nodes[1]),
                                        &nodes.at(element.nodes[2])};
    std::array<std::int64_t, 3> tags = {element.nodes[0], element.nodes[1], element.nodes[2]};

    const double signed_double_area = v[0]->x * (v[1]->y - v[2]->y) +
                                      v[1]->x * (v[2]->y - v[0]->y) + v[2]->x * (v[0]->y - v[1]->y);
    if (signed_double_area == 0.0)
      throw std::runtime_error("triangle mesh: zero-area element " + std::to_string(element.tag));
    if (signed_double_area < 0.0) {
      std::swap(v[1], v[2]);
      std::swap(tags[1], tags[2]);
    }

    elements.push_back(TriangleElement{
        tags[0], tags[1], tags[2], v[0]->x, v[0]->y, v[1]->x, v[1]->y, v[2]->x, v[2]->y,
        document.first_physical_tag(element.entity_dim, element.entity_tag)});
  }
  if (elements.empty())
    throw std::runtime_error("triangle mesh: document contains no triangle elements");
  return elements;
}

EdgeTable edge_adjacency(const std::vector<TriangleElement> &elements) {
  EdgeTable edges;
  edges.reserve(elements.size() * 3);
  for (std::size_t element = 0; element < elements.size(); ++element) {
    for (int local_edge = 0; local_edge < 3; ++local_edge) {
      const EdgeSide side{static_cast<std::int32_t>(element), local_edge,
                          elements[element].node(triangle_edge_vertices[local_edge][0]),
                          elements[element].node(triangle_edge_vertices[local_edge][1])};
      const FaceKey key = make_face_key(side.vertex0, side.vertex1);
      auto &slot = edges[key];
      if (slot.first.element == invalid_element) {
        slot.key = key;
        slot.first = side;
      } else if (slot.second.element == invalid_element) {
        slot.second = side;
      } else {
        throw std::runtime_error("triangle mesh: edge " + std::to_string(key.first) + "-" +
                                 std::to_string(key.second) +
                                 " has more than two adjacent elements (non-conforming mesh)");
      }
    }
  }
  return edges;
}

// Interior faces put the smaller element index on the left; boundary faces have
// only a left side. The normal comes from the left element's local edge.
TriangleFace make_face(const EdgeAdjacency &edge, const std::vector<TriangleElement> &elements,
                       const BoundaryEdgeTags &boundary_tags) {
  const bool interior = edge.second.element != invalid_element;
  const bool swap_sides = interior && edge.second.element < edge.first.element;
  const EdgeSide &left = swap_sides ? edge.second : edge.first;
  const EdgeSide right = interior ? (swap_sides ? edge.first : edge.second) : EdgeSide{};

  TriangleFace face{};
  face.vertex_a = edge.key.first;
  face.vertex_b = edge.key.second;
  face.left_element = left.element;
  face.right_element = right.element;
  face.left_local_edge = left.local_edge;
  face.right_local_edge = right.local_edge;
  face.left_vertex0 = left.vertex0;
  face.left_vertex1 = left.vertex1;
  face.right_vertex0 = right.vertex0;
  face.right_vertex1 = right.vertex1;
  elements[static_cast<std::size_t>(left.element)].map().physical_edge_normal(
      left.local_edge, face.outward_normal_x, face.outward_normal_y);
  face.physical_tag = interior ? invalid_physical : boundary_tags.at(edge.key);
  return face;
}

} // namespace

TriangleMesh build_triangle_mesh(const MshDocument &document) {
  const MshNodeLookup nodes(document, "triangle mesh");
  TriangleMesh mesh;
  mesh.physical_names = physical_names_of(document);
  mesh.elements = oriented_elements(document, nodes);
  for (std::size_t element = 0; element < mesh.elements.size(); ++element)
    mesh.physical_to_elements[mesh.elements[element].physical_tag].push_back(
        static_cast<std::int32_t>(element));

  const auto edges = edge_adjacency(mesh.elements);
  const BoundaryEdgeTags boundary_tags(document);
  mesh.faces.reserve(edges.size());
  for (const auto &[key, edge] : edges)
    mesh.faces.push_back(make_face(edge, mesh.elements, boundary_tags));

  // Hash-map iteration order is implementation-defined; sorting makes face
  // indices, flux summation order, and prepared-mesh output deterministic.
  std::sort(mesh.faces.begin(), mesh.faces.end(), [](const TriangleFace &a, const TriangleFace &b) {
    return a.left_element != b.left_element ? a.left_element < b.left_element
                                            : a.left_local_edge < b.left_local_edge;
  });

  for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
    const int tag = mesh.faces[face].physical_tag;
    if (tag != invalid_physical)
      mesh.physical_to_faces[tag].push_back(static_cast<std::int32_t>(face));
  }
  return mesh;
}

} // namespace wasl
