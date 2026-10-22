#include "wasl/mesh/prepared/prepared_mesh_validation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "wasl/mesh/mesh_index.hpp"

namespace wasl {
namespace {

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

bool all_finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

bool all_node_indices_valid(const std::vector<std::int64_t> &indices, std::size_t node_count) {
  return std::all_of(indices.begin(), indices.end(), [&](std::int64_t node) {
    return node >= 0 && static_cast<std::size_t>(node) < node_count;
  });
}

void validate_nodes(const PreparedMesh &mesh) {
  require(mesh.dimension >= 1 && mesh.dimension <= 2, "prepared mesh: dimension must be 1 or 2");
  require(mesh.nodes_per_element == mesh.dimension + 1,
          "prepared mesh: nodes_per_element does not match dimension");
  require(mesh.points.size() == mesh.node_tags.size() * 3,
          "prepared mesh: points size does not match node tags");
  require(all_finite(mesh.points), "prepared mesh: point coordinates must be finite");
}

void validate_elements(const PreparedMesh &mesh) {
  require(!mesh.element_nodes.empty() &&
              mesh.element_nodes.size() % static_cast<std::size_t>(mesh.nodes_per_element) == 0,
          "prepared mesh: invalid element connectivity");
  const auto element_count = mesh.element_nodes.size() / mesh.nodes_per_element;
  require(all_node_indices_valid(mesh.element_nodes, mesh.node_tags.size()),
          "prepared mesh: element has invalid local node index");
  require(mesh.element_physical_tags.size() == element_count,
          "prepared mesh: element physical-tag count mismatch");
  require(mesh.element_metrics.size() == element_count * mesh.metrics_per_element(),
          "prepared mesh: element metric count mismatch");
  require(all_finite(mesh.element_metrics), "prepared mesh: element metrics must be finite");
}

void validate_faces(const PreparedMesh &mesh) {
  const auto element_count = mesh.element_count();
  const auto face_count = mesh.face_count();
  require(mesh.face_nodes.size() == face_count * 2 && mesh.face_elements.size() == face_count * 2 &&
              mesh.face_local_edges.size() == face_count * 2 &&
              mesh.face_normals.size() == face_count * static_cast<std::size_t>(mesh.dimension) &&
              mesh.face_measures.size() == face_count,
          "prepared mesh: face array size mismatch");
  for (std::size_t face = 0; face < face_count; ++face) {
    const auto left = mesh.face_elements[face * 2];
    const auto right = mesh.face_elements[face * 2 + 1];
    require(left >= 0 && static_cast<std::size_t>(left) < element_count &&
                (right == invalid_element ||
                 (right >= 0 && static_cast<std::size_t>(right) < element_count)),
            "prepared mesh: face has invalid element index");
  }
  require(all_node_indices_valid(mesh.face_nodes, mesh.node_tags.size()),
          "prepared mesh: face has invalid local node index");
  require(all_finite(mesh.face_normals) &&
              std::all_of(mesh.face_measures.begin(), mesh.face_measures.end(),
                          [](double value) { return std::isfinite(value) && value > 0.0; }),
          "prepared mesh: face geometry must be finite with positive measures");
}

void validate_physical_names(const PreparedMesh &mesh) {
  require(!mesh.physical_name_offsets.empty() &&
              mesh.physical_name_dimensions.size() == mesh.physical_name_tags.size() &&
              mesh.physical_name_offsets.size() == mesh.physical_name_tags.size() + 1 &&
              mesh.physical_name_offsets.front() == 0 &&
              mesh.physical_name_offsets.back() == mesh.physical_name_bytes.size() &&
              std::is_sorted(mesh.physical_name_offsets.begin(), mesh.physical_name_offsets.end()),
          "prepared mesh: physical-name array size mismatch");
}

void validate_partition(const PreparedMesh &mesh) {
  const auto element_count = mesh.element_count();
  const auto rank_count = mesh.partition_rank_count;
  require(rank_count >= 1, "prepared mesh: partition_rank_count must be positive");
  require(mesh.element_owners.size() == element_count,
          "prepared mesh: element owner count mismatch");
  require(mesh.face_owners.size() == mesh.face_count(), "prepared mesh: face owner count mismatch");
  const auto in_range = [rank_count](std::int32_t owner) {
    return owner >= 0 && owner < rank_count;
  };
  require(std::all_of(mesh.element_owners.begin(), mesh.element_owners.end(), in_range),
          "prepared mesh: element owner is out of range");
  require(std::all_of(mesh.face_owners.begin(), mesh.face_owners.end(), in_range),
          "prepared mesh: face owner is out of range");

  // With at least as many elements as ranks, the contiguous block scheme gives
  // every rank work; with fewer, some ranks are legitimately idle.
  if (element_count >= static_cast<std::size_t>(rank_count)) {
    std::vector<unsigned char> owns_element(static_cast<std::size_t>(rank_count), 0);
    for (const auto owner : mesh.element_owners)
      owns_element[static_cast<std::size_t>(owner)] = 1;
    require(std::find(owns_element.begin(), owns_element.end(), 0) == owns_element.end(),
            "prepared mesh: a rank owns no elements");
  }

  for (std::size_t face = 0; face < mesh.face_count(); ++face) {
    const auto left = static_cast<std::size_t>(mesh.face_elements[face * 2]);
    require(mesh.face_owners[face] == mesh.element_owners.at(left),
            "prepared mesh: face owner does not match its left element owner");
  }
}

} // namespace

void validate_prepared_mesh(const PreparedMesh &mesh) {
  validate_nodes(mesh);
  validate_elements(mesh);
  validate_faces(mesh);
  validate_physical_names(mesh);
  validate_partition(mesh);
}

} // namespace wasl
