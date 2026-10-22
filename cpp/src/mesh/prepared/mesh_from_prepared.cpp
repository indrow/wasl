#include "wasl/mesh/prepared/mesh_from_prepared.hpp"

#include <cstddef>
#include <stdexcept>

namespace wasl {
namespace {

void require_single_rank(const PreparedMesh &prepared, int dimension, const char *message) {
  if (prepared.dimension != dimension || prepared.partition_rank_count != 1)
    throw std::runtime_error(message);
}

double point_x(const PreparedMesh &prepared, std::size_t node) { return prepared.points[node * 3]; }
double point_y(const PreparedMesh &prepared, std::size_t node) {
  return prepared.points[node * 3 + 1];
}

std::size_t element_node(const PreparedMesh &prepared, std::size_t element, int local_node) {
  return static_cast<std::size_t>(
      prepared.element_nodes[element * static_cast<std::size_t>(prepared.nodes_per_element) +
                             static_cast<std::size_t>(local_node)]);
}

} // namespace

LineMesh line_mesh_from_prepared(const PreparedMesh &prepared) {
  require_single_rank(prepared, 1, "simulation: CPU serial driver requires a single-rank 1D mesh");
  LineMesh mesh;
  for (std::size_t element = 0; element < prepared.element_count(); ++element) {
    const auto left = element_node(prepared, element, 0);
    const auto right = element_node(prepared, element, 1);
    mesh.elements.push_back({prepared.node_tags[left], prepared.node_tags[right],
                             point_x(prepared, left), point_x(prepared, right),
                             prepared.element_physical_tags[element]});
  }
  for (std::size_t face = 0; face < prepared.face_count(); ++face) {
    const auto node = static_cast<std::size_t>(prepared.face_nodes[face * 2]);
    mesh.faces.push_back({prepared.face_elements[face * 2], prepared.face_elements[face * 2 + 1],
                          point_x(prepared, node), prepared.face_normals[face],
                          prepared.face_physical_tags[face]});
  }
  return mesh;
}

TriangleMesh triangle_mesh_from_prepared(const PreparedMesh &prepared) {
  require_single_rank(prepared, 2, "simulation: CPU serial driver requires a single-rank 2D mesh");
  TriangleMesh mesh;
  for (std::size_t element = 0; element < prepared.element_count(); ++element) {
    const auto node0 = element_node(prepared, element, 0);
    const auto node1 = element_node(prepared, element, 1);
    const auto node2 = element_node(prepared, element, 2);
    mesh.elements.push_back({prepared.node_tags[node0], prepared.node_tags[node1],
                             prepared.node_tags[node2], point_x(prepared, node0),
                             point_y(prepared, node0), point_x(prepared, node1),
                             point_y(prepared, node1), point_x(prepared, node2),
                             point_y(prepared, node2), prepared.element_physical_tags[element]});
  }
  for (std::size_t face = 0; face < prepared.face_count(); ++face) {
    const auto left = prepared.face_elements[face * 2];
    const auto right = prepared.face_elements[face * 2 + 1];
    const auto left_edge = prepared.face_local_edges[face * 2];
    const auto right_edge = prepared.face_local_edges[face * 2 + 1];
    const auto &left_element = mesh.elements.at(static_cast<std::size_t>(left));
    TriangleFace triangle_face{
        prepared.node_tags[static_cast<std::size_t>(prepared.face_nodes[face * 2])],
        prepared.node_tags[static_cast<std::size_t>(prepared.face_nodes[face * 2 + 1])],
        left,
        right,
        left_edge,
        right_edge,
        left_element.node(triangle_edge_vertices[left_edge][0]),
        left_element.node(triangle_edge_vertices[left_edge][1]),
        0,
        0,
        prepared.face_normals[face * 2],
        prepared.face_normals[face * 2 + 1],
        prepared.face_physical_tags[face]};
    if (right != invalid_element) {
      const auto &right_element = mesh.elements.at(static_cast<std::size_t>(right));
      triangle_face.right_vertex0 = right_element.node(triangle_edge_vertices[right_edge][0]);
      triangle_face.right_vertex1 = right_element.node(triangle_edge_vertices[right_edge][1]);
    }
    mesh.faces.push_back(triangle_face);
  }
  return mesh;
}

} // namespace wasl
