// Deterministic partition ownership and one face-adjacent ghost-element layer.
//
// One MPI rank owns one partition of elements (ARCHITECTURE.md). Ownership and
// ghost sets are a function of prepared-mesh topology alone, so they are
// deterministic and testable without MPI.
//
// Element e is assigned to rank floor(e * rank_count / N), a contiguous block
// partition by element index. A face is owned by its left element's owner, the
// canonical face owner of NUMERICAL_CONVENTIONS.md.
//
// The DG residual and the WENO candidate transform touch only face-adjacent
// elements, so one ghost layer is sufficient: a rank imports every off-rank
// element that shares a face with one of its owned elements.

#pragma once

#include <cstdint>
#include <vector>

#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace wasl {

// Invariants:
//   rank_count >= 1.
//   element_owners and face_owners have one entry per element and face, each
//   in [0, rank_count).
//   Every rank in [0, min(rank_count, element_count)) owns an element.
//   face_owners[f] == element_owners[left element of f].
struct Partition {
  int rank_count = 1;
  std::vector<std::int32_t> element_owners;
  std::vector<std::int32_t> face_owners;
};

// Throws std::runtime_error if rank_count < 1.
Partition partition_mesh(const PreparedMesh &mesh, int rank_count);

// Record the partition's ownership arrays in the prepared mesh.
void assign_partition(PreparedMesh &mesh, const Partition &partition);

// Global indices of the elements a rank owns, ascending.
std::vector<std::int32_t> rank_owned_elements(const Partition &partition, int rank);

// Global indices of the off-rank, face-adjacent elements a rank must import,
// sorted and deduplicated. Boundary faces never contribute.
std::vector<std::int32_t> rank_ghost_elements(const PreparedMesh &mesh, const Partition &partition,
                                              int rank);

} // namespace wasl
