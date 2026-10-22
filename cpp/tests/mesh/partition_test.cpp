// Tests for deterministic partition ownership and one ghost-element layer.
//
// The fixture is a chain of six elements joined end to end through five
// interior faces, plus two boundary faces. Element and face indices are
// explicit so the partition math and ghost layer are unambiguous:
//
//   elements:  0 - 1 - 2 - 3 - 4 - 5
//   interior faces connect (i, i+1); boundary faces cap element 0 and 5.
//
// PreparedMesh stores faces as [left, right] with right == invalid_element on
// boundaries, so partition_mesh and rank_ghost_elements operate purely on
// element indices and face_elements. No geometry is needed for ownership.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_mesh.hpp" // invalid_element
#include "wasl/mesh/partition.hpp"

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

// A chain of element_count elements with element_count - 1 interior faces
// (each joining element i to i + 1) and two boundary faces. Face left/right
// follows the canonical convention: left is the smaller index.
wasl::PreparedMesh make_chain(std::size_t element_count) {
  wasl::PreparedMesh mesh;
  mesh.dimension = 1;
  mesh.nodes_per_element = 2;
  mesh.element_physical_tags.assign(element_count, 1);
  mesh.element_metrics.assign(element_count, 1.0);
  mesh.face_physical_tags.push_back(3); // left boundary of element 0
  mesh.face_elements.push_back(0);
  mesh.face_elements.push_back(wasl::invalid_element);
  mesh.face_nodes.insert(mesh.face_nodes.end(), {0, 0});
  mesh.face_local_edges.insert(mesh.face_local_edges.end(), {0, 0});
  mesh.face_normals.push_back(-1.0);
  mesh.face_measures.push_back(1.0);
  for (std::size_t i = 0; i + 1 < element_count; ++i) {
    mesh.face_physical_tags.push_back(wasl::invalid_physical);
    mesh.face_elements.push_back(static_cast<std::int32_t>(i));
    mesh.face_elements.push_back(static_cast<std::int32_t>(i + 1));
    mesh.face_nodes.insert(mesh.face_nodes.end(), {0, 0});
    mesh.face_local_edges.insert(mesh.face_local_edges.end(), {1, 0});
    mesh.face_normals.push_back(1.0);
    mesh.face_measures.push_back(1.0);
  }
  mesh.face_physical_tags.push_back(4); // right boundary of last element
  mesh.face_elements.push_back(static_cast<std::int32_t>(element_count - 1));
  mesh.face_elements.push_back(wasl::invalid_element);
  mesh.face_nodes.insert(mesh.face_nodes.end(), {0, 0});
  mesh.face_local_edges.insert(mesh.face_local_edges.end(), {1, 0});
  mesh.face_normals.push_back(1.0);
  mesh.face_measures.push_back(1.0);
  // Minimal physical-name bookkeeping so validate_prepared_mesh accepts it.
  mesh.physical_name_dimensions = {1, 1};
  mesh.physical_name_tags = {1, 3};
  mesh.physical_name_offsets = {0, 0, 0};
  mesh.node_tags = {1};
  mesh.points = {0.0, 0.0, 0.0};
  mesh.element_nodes.assign(element_count * 2, 0);
  // Backfill owners so the fixture is internally valid before partition tests
  // overwrite them.
  const auto partition = wasl::partition_mesh(mesh, 1);
  mesh.partition_rank_count = partition.rank_count;
  mesh.element_owners = partition.element_owners;
  mesh.face_owners = partition.face_owners;
  return mesh;
}

bool check_single_rank_is_whole_mesh() {
  auto mesh = make_chain(6);
  const auto partition = wasl::partition_mesh(mesh, 1);
  require(partition.rank_count == 1, "single rank: rank_count");
  require(partition.element_owners.size() == 6, "single rank: element owner count");
  require(partition.face_owners.size() == mesh.face_physical_tags.size(),
          "single rank: face owner count");
  for (const auto owner : partition.element_owners) {
    if (owner != 0)
      return false;
  }
  for (const auto owner : partition.face_owners) {
    if (owner != 0)
      return false;
  }
  // One rank owns everything, so it imports no ghosts.
  const auto ghosts = wasl::rank_ghost_elements(mesh, partition, 0);
  return ghosts.empty();
}

bool check_two_rank_block_and_face_ownership() {
  const auto mesh = make_chain(6);
  const auto partition = wasl::partition_mesh(mesh, 2);
  // floor(e * 2 / 6): {0,0,0,1,1,1}.
  const std::vector<std::int32_t> expected_owners = {0, 0, 0, 1, 1, 1};
  if (partition.element_owners != expected_owners)
    return false;
  // Every face owner equals its left element's owner.
  for (std::size_t f = 0; f < mesh.face_physical_tags.size(); ++f) {
    const auto left = mesh.face_elements[f * 2];
    if (partition.face_owners[f] != partition.element_owners[static_cast<std::size_t>(left)])
      return false;
  }
  return true;
}

bool check_ghost_layer_is_face_adjacent_only() {
  const auto mesh = make_chain(6);
  const auto partition = wasl::partition_mesh(mesh, 2);
  // Rank 0 owns {0,1,2}. The only cross-rank interior face is (2,3), so it
  // imports exactly element 3. Element 4 is a neighbor's neighbor (second ring).
  const auto ghosts0 = wasl::rank_ghost_elements(mesh, partition, 0);
  if (ghosts0 != std::vector<std::int32_t>{3})
    return false;
  // Rank 1 owns {3,4,5}. Its only cross-rank interior face is (2,3); it imports
  // exactly element 2.
  const auto ghosts1 = wasl::rank_ghost_elements(mesh, partition, 1);
  if (ghosts1 != std::vector<std::int32_t>{2})
    return false;
  // Ghosts are always sorted and deduplicated.
  return std::is_sorted(ghosts0.begin(), ghosts0.end()) &&
         std::is_sorted(ghosts1.begin(), ghosts1.end());
}

// Acceptance check: the two rank-local slices plus ghost layers reconstruct the
// serial face topology with no duplicate owned elements. Concretely, across all
// ranks the owned element sets are disjoint and cover every element, and every
// interior face reappears once as an owned side and once as a ghost side.
bool check_two_rank_reconstructs_serial_topology() {
  const auto mesh = make_chain(6);
  const auto partition = wasl::partition_mesh(mesh, 2);

  // Every element is owned by exactly one rank, and ranks cover all elements.
  std::vector<std::int32_t> owned_union(6, -1);
  for (std::size_t e = 0; e < 6; ++e) {
    owned_union[e] = partition.element_owners[e];
  }
  std::set<std::int32_t> owner_set(owned_union.begin(), owned_union.end());
  if (owner_set.count(0) == 0 || owner_set.count(1) == 0)
    return false;

  // For every interior face, one endpoint is owned and the other is a ghost on
  // the owning rank. Reconstruct each face as an (owned, ghost) pair and show
  // they match the serial (left, right) pairs.
  for (std::size_t f = 0; f < mesh.face_physical_tags.size(); ++f) {
    const auto left = mesh.face_elements[f * 2];
    const auto right = mesh.face_elements[f * 2 + 1];
    if (right == wasl::invalid_element)
      continue; // boundary face has no neighbor to reconstruct
    const auto left_rank = partition.element_owners[static_cast<std::size_t>(left)];
    const auto right_rank = partition.element_owners[static_cast<std::size_t>(right)];
    if (left_rank == right_rank)
      continue; // same-rank interior face is fully local, no ghost needed
    const auto ghosts = wasl::rank_ghost_elements(mesh, partition, left_rank);
    if (std::find(ghosts.begin(), ghosts.end(), right) == ghosts.end())
      return false;
    const auto other_ghosts = wasl::rank_ghost_elements(mesh, partition, right_rank);
    if (std::find(other_ghosts.begin(), other_ghosts.end(), left) == other_ghosts.end())
      return false;
  }
  // No element is owned by more than one rank (already implied by element_owners
  // being a function), and the two owned sets are disjoint.
  std::vector<std::int32_t> rank0;
  std::vector<std::int32_t> rank1;
  for (std::size_t e = 0; e < 6; ++e) {
    (partition.element_owners[e] == 0 ? rank0 : rank1).push_back(static_cast<std::int32_t>(e));
  }
  std::set<std::int32_t> overlap;
  std::set_intersection(rank0.begin(), rank0.end(), rank1.begin(), rank1.end(),
                        std::inserter(overlap, overlap.end()));
  return overlap.empty();
}

bool check_three_rank_balanced_split() {
  const auto mesh = make_chain(6);
  const auto partition = wasl::partition_mesh(mesh, 3);
  // floor(e * 3 / 6): {0,0,1,1,2,2}.
  if (partition.element_owners != std::vector<std::int32_t>{0, 0, 1, 1, 2, 2})
    return false;
  // Rank 1 owns {2,3}; it straddles two partition boundaries, (1,2) and (3,4),
  // importing elements 1 and 4.
  const auto ghosts = wasl::rank_ghost_elements(mesh, partition, 1);
  return ghosts == std::vector<std::int32_t>{1, 4};
}

bool check_rejects_nonpositive_rank_count() {
  auto mesh = make_chain(4);
  try {
    static_cast<void>(wasl::partition_mesh(mesh, 0));
  } catch (const std::runtime_error &error) {
    return std::string{error.what()}.find("rank_count") != std::string::npos;
  }
  return false;
}

} // namespace

int main() {
  try {
    if (!check_single_rank_is_whole_mesh())
      return 1;
    if (!check_two_rank_block_and_face_ownership())
      return 1;
    if (!check_ghost_layer_is_face_adjacent_only())
      return 1;
    if (!check_two_rank_reconstructs_serial_topology())
      return 1;
    if (!check_three_rank_balanced_split())
      return 1;
    if (!check_rejects_nonpositive_rank_count())
      return 1;
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
