#include "wasl/mesh/prepared/prepare_mesh.hpp"

#include <cstddef>

#include "mesh/prepared/prepared_node_table.hpp"
#include "wasl/mesh/partition/partition.hpp"

namespace wasl {
namespace {

void pack_physical_names(PreparedMesh &result, const PhysicalNames &names) {
  result.physical_name_offsets.push_back(0);
  for (const auto &[key, name] : names) {
    result.physical_name_dimensions.push_back(key.first);
    result.physical_name_tags.push_back(key.second);
    result.physical_name_bytes.insert(result.physical_name_bytes.end(), name.begin(), name.end());
    result.physical_name_offsets.push_back(result.physical_name_bytes.size());
  }
}

// Multi-rank partitioning is applied later, once a rank count is known
// (ARCHITECTURE.md).
void assign_single_rank(PreparedMesh &result) {
  assign_partition(result, partition_mesh(result, 1));
}

} // namespace

PreparedMesh prepare_mesh(const LineMesh &mesh) {
  PreparedMesh result;
  result.dimension = 1;
  result.nodes_per_element = 2;
  PreparedNodeTable nodes(result);
  for (const auto &element : mesh.elements) {
    const auto left = nodes.add(element.left_node, element.x_left, 0.0);
    const auto right = nodes.add(element.right_node, element.x_right, 0.0);
    result.element_nodes.insert(result.element_nodes.end(), {left, right});
    result.element_physical_tags.push_back(element.physical_tag);
    result.element_metrics.push_back((element.x_right - element.x_left) / 2.0);
  }
  for (const auto &face : mesh.faces) {
    const auto &element = mesh.elements.at(static_cast<std::size_t>(face.left_element));
    const auto node = nodes.at(face.x == element.x_left ? element.left_node : element.right_node);
    result.face_nodes.insert(result.face_nodes.end(), {node, node});
    result.face_elements.insert(result.face_elements.end(),
                                {face.left_element, face.right_element});
    result.face_local_edges.insert(result.face_local_edges.end(),
                                   {face.outward_normal < 0.0 ? 0 : 1, 0});
    result.face_normals.push_back(face.outward_normal);
    result.face_measures.push_back(1.0);
    result.face_physical_tags.push_back(face.physical_tag);
  }
  pack_physical_names(result, mesh.physical_names);
  assign_single_rank(result);
  return result;
}

PreparedMesh prepare_mesh(const TriangleMesh &mesh) {
  PreparedMesh result;
  result.dimension = 2;
  result.nodes_per_element = 3;
  PreparedNodeTable nodes(result);
  for (const auto &element : mesh.elements) {
    const auto node0 = nodes.add(element.node0, element.x0, element.y0);
    const auto node1 = nodes.add(element.node1, element.x1, element.y1);
    const auto node2 = nodes.add(element.node2, element.x2, element.y2);
    result.element_nodes.insert(result.element_nodes.end(), {node0, node1, node2});
    result.element_physical_tags.push_back(element.physical_tag);
    const auto map = element.map();
    result.element_metrics.insert(result.element_metrics.end(),
                                  {map.ax() / 2.0, map.bx() / 2.0, map.ay() / 2.0, map.by() / 2.0,
                                   map.jacobian_determinant()});
  }
  for (const auto &face : mesh.faces) {
    result.face_nodes.insert(result.face_nodes.end(),
                             {nodes.at(face.vertex_a), nodes.at(face.vertex_b)});
    result.face_elements.insert(result.face_elements.end(),
                                {face.left_element, face.right_element});
    result.face_local_edges.insert(result.face_local_edges.end(),
                                   {face.left_local_edge, face.right_local_edge});
    result.face_normals.insert(result.face_normals.end(),
                               {face.outward_normal_x, face.outward_normal_y});
    const auto &owner = mesh.elements.at(static_cast<std::size_t>(face.left_element));
    result.face_measures.push_back(owner.map().edge_length(face.left_local_edge));
    result.face_physical_tags.push_back(face.physical_tag);
  }
  pack_physical_names(result, mesh.physical_names);
  assign_single_rank(result);
  return result;
}

} // namespace wasl
