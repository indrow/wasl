#include "wasl/case/prepared_mesh.hpp"

#include <stdexcept>

namespace wasl {

LineMesh line_mesh_from_prepared(const PreparedMesh &prepared) {
  if (prepared.dimension != 1 || prepared.partition_rank_count != 1) {
    throw std::runtime_error("simulation: CPU serial driver requires a single-rank 1D mesh");
  }
  LineMesh mesh;
  for (std::size_t element = 0; element < prepared.element_physical_tags.size(); ++element) {
    const auto left = static_cast<std::size_t>(prepared.element_nodes[element * 2]);
    const auto right = static_cast<std::size_t>(prepared.element_nodes[element * 2 + 1]);
    mesh.elements.push_back({prepared.node_tags[left], prepared.node_tags[right],
                             prepared.points[left * 3], prepared.points[right * 3],
                             prepared.element_physical_tags[element]});
  }
  for (std::size_t face = 0; face < prepared.face_physical_tags.size(); ++face) {
    const auto node = static_cast<std::size_t>(prepared.face_nodes[face * 2]);
    mesh.faces.push_back({prepared.face_elements[face * 2], prepared.face_elements[face * 2 + 1],
                          prepared.points[node * 3], prepared.face_normals[face],
                          prepared.face_physical_tags[face]});
  }
  return mesh;
}

TriangleMesh triangle_mesh_from_prepared(const PreparedMesh &prepared) {
  if (prepared.dimension != 2 || prepared.partition_rank_count != 1) {
    throw std::runtime_error("simulation: CPU serial driver requires a single-rank 2D mesh");
  }
  TriangleMesh mesh;
  for (std::size_t element = 0; element < prepared.element_physical_tags.size(); ++element) {
    const auto node0 = static_cast<std::size_t>(prepared.element_nodes[element * 3]);
    const auto node1 = static_cast<std::size_t>(prepared.element_nodes[element * 3 + 1]);
    const auto node2 = static_cast<std::size_t>(prepared.element_nodes[element * 3 + 2]);
    mesh.elements.push_back(
        {prepared.node_tags[node0], prepared.node_tags[node1], prepared.node_tags[node2],
         prepared.points[node0 * 3], prepared.points[node0 * 3 + 1], prepared.points[node1 * 3],
         prepared.points[node1 * 3 + 1], prepared.points[node2 * 3], prepared.points[node2 * 3 + 1],
         prepared.element_physical_tags[element]});
  }
  constexpr int edge_vertices[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  for (std::size_t face = 0; face < prepared.face_physical_tags.size(); ++face) {
    const auto left = prepared.face_elements[face * 2];
    const auto right = prepared.face_elements[face * 2 + 1];
    const auto left_edge = prepared.face_local_edges[face * 2];
    const auto right_edge = prepared.face_local_edges[face * 2 + 1];
    const auto node_a = static_cast<std::size_t>(prepared.face_nodes[face * 2]);
    const auto node_b = static_cast<std::size_t>(prepared.face_nodes[face * 2 + 1]);
    const auto &left_element = mesh.elements.at(static_cast<std::size_t>(left));
    TriangleFace triangle_face{prepared.node_tags[node_a],
                               prepared.node_tags[node_b],
                               left,
                               right,
                               left_edge,
                               right_edge,
                               left_element.node(edge_vertices[left_edge][0]),
                               left_element.node(edge_vertices[left_edge][1]),
                               0,
                               0,
                               prepared.face_normals[face * 2],
                               prepared.face_normals[face * 2 + 1],
                               prepared.face_physical_tags[face]};
    if (right != invalid_element) {
      const auto &right_element = mesh.elements.at(static_cast<std::size_t>(right));
      triangle_face.right_vertex0 = right_element.node(edge_vertices[right_edge][0]);
      triangle_face.right_vertex1 = right_element.node(edge_vertices[right_edge][1]);
    }
    mesh.faces.push_back(triangle_face);
  }
  return mesh;
}

} // namespace wasl
