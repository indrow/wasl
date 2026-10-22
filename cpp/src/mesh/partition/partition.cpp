#include "wasl/mesh/partition/partition.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

#include "wasl/mesh/mesh_index.hpp"

namespace wasl {

Partition partition_mesh(const PreparedMesh &mesh, int rank_count) {
  if (rank_count < 1)
    throw std::runtime_error("partition: rank_count must be positive");
  const auto element_count = mesh.element_count();

  Partition partition;
  partition.rank_count = rank_count;
  partition.element_owners.reserve(element_count);
  // 64-bit products avoid overflow on large meshes.
  for (std::size_t element = 0; element < element_count; ++element)
    partition.element_owners.push_back(static_cast<std::int32_t>(
        (static_cast<std::uint64_t>(element) * static_cast<std::uint64_t>(rank_count)) /
        element_count));

  partition.face_owners.reserve(mesh.face_count());
  for (std::size_t face = 0; face < mesh.face_count(); ++face) {
    const auto left = static_cast<std::size_t>(mesh.face_elements[face * 2]);
    partition.face_owners.push_back(partition.element_owners.at(left));
  }
  return partition;
}

void assign_partition(PreparedMesh &mesh, const Partition &partition) {
  mesh.partition_rank_count = partition.rank_count;
  mesh.element_owners = partition.element_owners;
  mesh.face_owners = partition.face_owners;
}

std::vector<std::int32_t> rank_owned_elements(const Partition &partition, int rank) {
  std::vector<std::int32_t> owned;
  for (std::size_t element = 0; element < partition.element_owners.size(); ++element) {
    if (partition.element_owners[element] == rank)
      owned.push_back(static_cast<std::int32_t>(element));
  }
  return owned;
}

std::vector<std::int32_t> rank_ghost_elements(const PreparedMesh &mesh, const Partition &partition,
                                              int rank) {
  std::vector<std::int32_t> ghosts;
  std::vector<unsigned char> is_ghost(mesh.element_count(), 0);
  const auto add_ghost = [&](std::int32_t element) {
    auto &flag = is_ghost[static_cast<std::size_t>(element)];
    if (!flag) {
      flag = 1;
      ghosts.push_back(element);
    }
  };

  for (std::size_t face = 0; face < mesh.face_count(); ++face) {
    const auto left = mesh.face_elements[face * 2];
    const auto right = mesh.face_elements[face * 2 + 1];
    if (right == invalid_element)
      continue;
    const auto left_rank = partition.element_owners.at(static_cast<std::size_t>(left));
    const auto right_rank = partition.element_owners.at(static_cast<std::size_t>(right));
    if (left_rank == rank && right_rank != rank)
      add_ghost(right);
    else if (right_rank == rank && left_rank != rank)
      add_ghost(left);
  }
  std::sort(ghosts.begin(), ghosts.end());
  return ghosts;
}

} // namespace wasl
