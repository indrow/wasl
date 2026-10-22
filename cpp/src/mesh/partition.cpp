#include "wasl/mesh/partition.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp" // invalid_element

namespace wasl {

Partition partition_mesh(const PreparedMesh &mesh, int rank_count) {
  if (rank_count < 1) {
    throw std::runtime_error("partition: rank_count must be positive");
  }
  const auto element_count = mesh.element_physical_tags.size();
  const auto face_count = mesh.face_physical_tags.size();

  Partition partition;
  partition.rank_count = rank_count;
  partition.element_owners.reserve(element_count);
  // Contiguous block: element e -> rank floor(e * rank_count / N). Iterating
  // with 64-bit products avoids overflow on large meshes.
  for (std::size_t e = 0; e < element_count; ++e) {
    const auto rank = static_cast<std::int32_t>(
        (static_cast<std::uint64_t>(e) * static_cast<std::uint64_t>(rank_count)) / element_count);
    partition.element_owners.push_back(rank);
  }

  // A face is owned by its left element's owner. face_elements is stored as
  // [left, right] per face with right == invalid_element on boundaries.
  partition.face_owners.reserve(face_count);
  for (std::size_t f = 0; f < face_count; ++f) {
    const auto left = mesh.face_elements[f * 2];
    partition.face_owners.push_back(partition.element_owners.at(static_cast<std::size_t>(left)));
  }
  return partition;
}

std::vector<std::int32_t> rank_ghost_elements(const PreparedMesh &mesh, const Partition &partition,
                                              int rank) {
  const auto element_count = mesh.element_physical_tags.size();
  const auto face_count = mesh.face_physical_tags.size();
  std::vector<std::int32_t> ghosts;
  std::vector<unsigned char> is_ghost(element_count, 0);

  for (std::size_t f = 0; f < face_count; ++f) {
    const auto left = mesh.face_elements[f * 2];
    const auto right = mesh.face_elements[f * 2 + 1];
    if (right == invalid_element) {
      continue; // boundary face has no neighbor
    }
    const auto left_rank = partition.element_owners.at(static_cast<std::size_t>(left));
    const auto right_rank = partition.element_owners.at(static_cast<std::size_t>(right));
    // Only faces straddling a partition boundary contribute off-rank neighbors.
    if (left_rank == rank && right_rank != rank) {
      const auto g = static_cast<std::size_t>(right);
      if (!is_ghost[g]) {
        is_ghost[g] = 1;
        ghosts.push_back(right);
      }
    } else if (right_rank == rank && left_rank != rank) {
      const auto g = static_cast<std::size_t>(left);
      if (!is_ghost[g]) {
        is_ghost[g] = 1;
        ghosts.push_back(left);
      }
    }
  }
  std::sort(ghosts.begin(), ghosts.end());
  return ghosts;
}

} // namespace wasl
