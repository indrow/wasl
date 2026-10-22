#include "wasl/mesh/triangle_mesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "wasl/mesh/msh_reader.hpp"

namespace wasl {
namespace {

using NodeLookup = std::unordered_map<std::int64_t, const MshNode *>;

NodeLookup index_nodes(const MshDocument &document) {
  NodeLookup nodes;
  nodes.reserve(document.nodes.size());
  for (const auto &node : document.nodes) {
    if (!nodes.emplace(node.tag, &node).second) {
      throw std::runtime_error("triangle mesh: duplicate node tag " + std::to_string(node.tag));
    }
  }
  return nodes;
}

const MshNode &require_node(const NodeLookup &nodes, std::int64_t tag) {
  const auto found = nodes.find(tag);
  if (found == nodes.end()) {
    throw std::runtime_error("triangle mesh: element references missing node " +
                             std::to_string(tag));
  }
  return *found->second;
}

// The three local edges of a triangle, as (start local vertex, end local vertex).
constexpr std::array<std::array<int, 2>, 3> kLocalEdges = {{{0, 1}, {1, 2}, {2, 0}}};

struct FaceAccumulator {
  std::int64_t key_lo = 0; // smaller vertex ID of the sorted key
  std::int64_t key_hi = 0; // larger vertex ID of the sorted key
  std::int32_t first_element = invalid_element;
  int first_local_edge = 0;
  std::int64_t first_v0 = 0; // ordered start vertex as first element traverses the edge
  std::int64_t first_v1 = 0; // ordered end vertex
  std::int32_t second_element = invalid_element; // invalid_element until a second adjacency
  int second_local_edge = 0;
  std::int64_t second_v0 = 0;
  std::int64_t second_v1 = 0;
};

FaceKey make_key(std::int64_t a, std::int64_t b) { return a < b ? FaceKey{a, b} : FaceKey{b, a}; }

// Resolve the boundary physical tag for a boundary edge whose sorted vertex key
// is (a, b). A boundary edge lies on one curve entity whose two bounding points
// match these vertices; the curve's first physical tag names the boundary. Gmsh
// curve bounding points may be signed (sign encodes orientation); only the
// magnitude is the node tag.
// Boundary line elements map each meshed edge (sorted node pair) to the
// physical tag of the curve entity that owns it. A curve split into many
// segments only has its end nodes in the entity bounding points, so the
// element map is authoritative whenever line elements are present.
using EdgeTags = std::unordered_map<FaceKey, int, FaceKeyHash>;

EdgeTags boundary_edge_tags(const MshDocument &document) {
  EdgeTags tags;
  for (const auto &element : document.elements) {
    if (element.element_type != 1 || element.nodes.size() != 2) {
      continue;
    }
    tags.emplace(make_key(element.nodes[0], element.nodes[1]),
                 document.first_physical_tag(element.entity_dim, element.entity_tag));
  }
  return tags;
}

int boundary_physical_tag(const MshDocument &document, const EdgeTags &edge_tags, std::int64_t a,
                          std::int64_t b) {
  const auto found = edge_tags.find(FaceKey{a, b});
  if (found != edge_tags.end()) {
    return found->second;
  }
  for (const auto &curve : document.curve_entities) {
    if (curve.bounding_points.size() != 2) {
      continue;
    }
    const auto p0 = static_cast<std::int64_t>(std::abs(curve.bounding_points[0]));
    const auto p1 = static_cast<std::int64_t>(std::abs(curve.bounding_points[1]));
    const FaceKey key = make_key(p0, p1);
    if (key.first == a && key.second == b) {
      return curve.physical_tags.empty() ? invalid_physical : curve.physical_tags.front();
    }
  }
  return invalid_physical;
}

} // namespace

TriangleMesh build_triangle_mesh(const MshDocument &document) {
  TriangleMesh mesh;
  for (const auto &physical : document.physical_names) {
    mesh.physical_names[{physical.dim, physical.tag}] = physical.name;
  }

  const auto nodes = index_nodes(document);

  // Collect and orient triangle elements. Vertex order is normalized to positive
  // (counterclockwise) physical orientation by swapping the trailing two
  // vertices when the signed area is negative.
  mesh.elements.reserve(document.elements.size());
  for (const auto &element : document.elements) {
    if (element.element_type != 2 || element.nodes.size() != 3) {
      continue;
    }
    std::array<const MshNode *, 3> v = {
        &require_node(nodes, element.nodes[0]),
        &require_node(nodes, element.nodes[1]),
        &require_node(nodes, element.nodes[2]),
    };
    std::array<std::int64_t, 3> tags = {element.nodes[0], element.nodes[1], element.nodes[2]};

    const double signed_double = v[0]->x * (v[1]->y - v[2]->y) + v[1]->x * (v[2]->y - v[0]->y) +
                                 v[2]->x * (v[0]->y - v[1]->y);
    if (signed_double == 0.0) {
      throw std::runtime_error("triangle mesh: zero-area element " + std::to_string(element.tag));
    }
    if (signed_double < 0.0) {
      std::swap(v[1], v[2]);
      std::swap(tags[1], tags[2]);
    }

    mesh.elements.push_back(TriangleElement{
        tags[0],
        tags[1],
        tags[2],
        v[0]->x,
        v[0]->y,
        v[1]->x,
        v[1]->y,
        v[2]->x,
        v[2]->y,
        document.first_physical_tag(element.entity_dim, element.entity_tag),
    });
  }
  if (mesh.elements.empty()) {
    throw std::runtime_error("triangle mesh: document contains no triangle elements");
  }

  for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
    mesh.physical_to_elements[mesh.elements[e].physical_tag].push_back(
        static_cast<std::int32_t>(e));
  }

  // Build edges. Key each by the sorted pair of global vertex IDs and accumulate
  // adjacencies. Each element contributes three local edges; record the ordered
  // vertex pair as the element traverses its local edge.
  std::unordered_map<FaceKey, FaceAccumulator, FaceKeyHash> edges;
  edges.reserve(mesh.elements.size() * 3);
  for (std::int32_t e = 0; e < static_cast<std::int32_t>(mesh.elements.size()); ++e) {
    const auto &elem = mesh.elements[static_cast<std::size_t>(e)];
    for (int le = 0; le < 3; ++le) {
      const int va = kLocalEdges[static_cast<std::size_t>(le)][0];
      const int vb = kLocalEdges[static_cast<std::size_t>(le)][1];
      const std::int64_t g0 = elem.node(va);
      const std::int64_t g1 = elem.node(vb);
      const FaceKey key = make_key(g0, g1);
      auto &slot = edges[key];
      if (slot.first_element == invalid_element) {
        slot.key_lo = key.first;
        slot.key_hi = key.second;
        slot.first_element = e;
        slot.first_local_edge = le;
        slot.first_v0 = g0;
        slot.first_v1 = g1;
      } else if (slot.second_element == invalid_element) {
        slot.second_element = e;
        slot.second_local_edge = le;
        slot.second_v0 = g0;
        slot.second_v1 = g1;
      } else {
        throw std::runtime_error("triangle mesh: edge " + std::to_string(key.first) + "-" +
                                 std::to_string(key.second) +
                                 " has more than two adjacent elements (non-conforming mesh)");
      }
    }
  }

  const auto edge_tags = boundary_edge_tags(document);

  // Emit faces with deterministic left/right ownership: for interior faces the
  // element with the smaller index is left, the other is right, and per-side
  // vertex order records each element's traversal. Boundary faces carry an
  // outward physical normal.
  mesh.faces.reserve(edges.size());
  for (auto &[key, acc] : edges) {
    const bool interior = acc.second_element != invalid_element;
    TriangleFace face{};
    face.vertex_a = acc.key_lo;
    face.vertex_b = acc.key_hi;
    if (interior) {
      const bool ordered = acc.first_element < acc.second_element;
      face.left_element = ordered ? acc.first_element : acc.second_element;
      face.right_element = ordered ? acc.second_element : acc.first_element;
      face.left_local_edge = ordered ? acc.first_local_edge : acc.second_local_edge;
      face.right_local_edge = ordered ? acc.second_local_edge : acc.first_local_edge;
      face.left_vertex0 = ordered ? acc.first_v0 : acc.second_v0;
      face.left_vertex1 = ordered ? acc.first_v1 : acc.second_v1;
      face.right_vertex0 = ordered ? acc.second_v0 : acc.first_v0;
      face.right_vertex1 = ordered ? acc.second_v1 : acc.first_v1;
      face.physical_tag = invalid_physical;
      // Canonical interior normal: from left element toward right element, taken
      // from the left element geometry along its owning local edge.
      const auto &owner = mesh.elements[static_cast<std::size_t>(face.left_element)];
      const auto map = TriangleMap{owner.x0, owner.y0, owner.x1, owner.y1, owner.x2, owner.y2};
      map.physical_edge_normal(face.left_local_edge, face.outward_normal_x, face.outward_normal_y);
    } else {
      face.left_element = acc.first_element;
      face.right_element = invalid_element;
      face.left_local_edge = acc.first_local_edge;
      face.right_local_edge = 0;
      face.left_vertex0 = acc.first_v0;
      face.left_vertex1 = acc.first_v1;
      face.right_vertex0 = 0;
      face.right_vertex1 = 0;
      const auto &owner = mesh.elements[static_cast<std::size_t>(face.left_element)];
      const auto map = TriangleMap{owner.x0, owner.y0, owner.x1, owner.y1, owner.x2, owner.y2};
      map.physical_edge_normal(face.left_local_edge, face.outward_normal_x, face.outward_normal_y);
      face.physical_tag = boundary_physical_tag(document, edge_tags, acc.key_lo, acc.key_hi);
    }
    mesh.faces.push_back(std::move(face));
  }
  // Hash-map iteration order is implementation-defined; sort so face indices,
  // flux summation order, and prepared-mesh output are deterministic.
  std::sort(mesh.faces.begin(), mesh.faces.end(),
            [](const TriangleFace &left, const TriangleFace &right) {
              return left.left_element != right.left_element
                         ? left.left_element < right.left_element
                         : left.left_local_edge < right.left_local_edge;
            });

  for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
    const int tag = mesh.faces[i].physical_tag;
    if (tag != invalid_physical) {
      mesh.physical_to_faces[tag].push_back(static_cast<std::int32_t>(i));
    }
  }
  return mesh;
}

TriangleMap triangle_element_map(const TriangleMesh &mesh, std::size_t element) {
  const auto &e = mesh.elements.at(element);
  return TriangleMap{e.x0, e.y0, e.x1, e.y1, e.x2, e.y2};
}

std::size_t count_interior_faces(const TriangleMesh &mesh) noexcept {
  return static_cast<std::size_t>(
      std::count_if(mesh.faces.begin(), mesh.faces.end(), [](const TriangleFace &face) {
        return face.right_element != invalid_element;
      }));
}

std::size_t count_boundary_faces(const TriangleMesh &mesh) noexcept {
  return mesh.faces.size() - count_interior_faces(mesh);
}

} // namespace wasl
