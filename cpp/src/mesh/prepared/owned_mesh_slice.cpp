#include "wasl/mesh/prepared/owned_mesh_slice.hpp"

#include <cstddef>
#include <stdexcept>

#include "mesh/prepared/prepared_node_table.hpp"
#include "wasl/mesh/prepared/prepared_mesh_validation.hpp"

namespace wasl {

PreparedMesh owned_mesh_slice(const PreparedMesh &mesh, int rank) {
  validate_prepared_mesh(mesh);
  if (rank < 0 || rank >= mesh.partition_rank_count)
    throw std::runtime_error("owned mesh slice: rank is out of range");

  PreparedMesh result;
  result.dimension = mesh.dimension;
  result.nodes_per_element = mesh.nodes_per_element;
  result.physical_name_dimensions = mesh.physical_name_dimensions;
  result.physical_name_tags = mesh.physical_name_tags;
  result.physical_name_offsets = mesh.physical_name_offsets;
  result.physical_name_bytes = mesh.physical_name_bytes;

  PreparedNodeTable nodes(result);
  const auto nodes_per_element = static_cast<std::size_t>(mesh.nodes_per_element);
  const auto metrics_per_element = mesh.metrics_per_element();
  for (std::size_t element = 0; element < mesh.element_owners.size(); ++element) {
    if (mesh.element_owners[element] != rank)
      continue;
    for (std::size_t local_node = 0; local_node < nodes_per_element; ++local_node) {
      const auto source =
          static_cast<std::size_t>(mesh.element_nodes[element * nodes_per_element + local_node]);
      result.element_nodes.push_back(
          nodes.add(mesh.node_tags[source], mesh.points[source * 3], mesh.points[source * 3 + 1]));
    }
    result.element_physical_tags.push_back(mesh.element_physical_tags[element]);
    result.element_metrics.insert(
        result.element_metrics.end(), mesh.element_metrics.begin() + element * metrics_per_element,
        mesh.element_metrics.begin() + (element + 1) * metrics_per_element);
    result.element_owners.push_back(0);
  }
  if (result.element_nodes.empty())
    throw std::runtime_error("owned mesh slice: rank owns no elements");
  return result;
}

} // namespace wasl
