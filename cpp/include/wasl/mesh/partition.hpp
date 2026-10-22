// Deterministic partition ownership and one face-adjacent ghost-element layer.
//
// One MPI rank owns one partition of elements (ARCHITECTURE.md). This header computes
// the metadata only: which rank owns each element and face, and which off-rank
// elements each rank must import. The transport (host-staged MPI) is a later
// task; ownership and ghost sets are a function of the prepared-mesh topology
// alone, so they are deterministic and unit-testable without MPI.
//
// Ownership scheme. Element e is assigned to rank floor(e * rank_count / N),
// a contiguous block partition by element index. It is reproducible, cheap,
// and free of external dependencies (ARCHITECTURE.md). A face is owned by its left
// element's owner, matching the canonical face-owner convention in
// NUMERICAL_CONVENTIONS.md: the left element is the deterministic face owner.
//
// Ghost layer. The DG residual and the WENO candidate transform touch only the
// elements directly adjacent to a face (left/right), never a second ring. So
// one face-adjacent layer is exactly sufficient: a rank imports every off-rank
// element that shares a face with one of its owned elements, deduplicated.

#pragma once

#include <cstdint>
#include <vector>

namespace wasl {

struct PreparedMesh; // defined in wasl/io/case_io.hpp

// Per-rank partition metadata over the whole mesh.
//
// Invariants:
//   rank_count >= 1.
//   element_owners and face_owners have one entry per prepared-mesh element and
//   face; every value is in [0, rank_count).
//   Every rank in [0, min(rank_count, element_count)) owns at least one element.
//   face_owners[f] == element_owners[face left element of f].
struct Partition {
  int rank_count = 1;
  std::vector<std::int32_t> element_owners; // [element] owning rank
  std::vector<std::int32_t> face_owners;    // [face] owning rank
};

// Assign each element to a contiguous block of ranks and derive face ownership.
// Throws std::runtime_error if rank_count < 1.
Partition partition_mesh(const PreparedMesh &mesh, int rank_count);

// Global indices of the off-rank elements one rank must import: the elements
// face-adjacent to its owned elements that another rank owns. The result is
// sorted and deduplicated. Boundary faces (no right element) never contribute.
std::vector<std::int32_t> rank_ghost_elements(const PreparedMesh &mesh, const Partition &partition,
                                              int rank);

} // namespace wasl
