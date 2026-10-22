// Flat, file-ready mesh representation shared by preparation, partitioning,
// HDF5 I/O, and visualization.
//
// Nodes use local indices into node_tags; points are packed (x, y, z). Element
// connectivity holds nodes_per_element local node indices per element. Faces
// store [left, right] element pairs (right == invalid_element on a boundary),
// [left, right] local edges, and one normal component per dimension. Physical
// names are packed bytes with offsets.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace wasl {

struct PreparedMesh {
  int dimension = 0;
  int nodes_per_element = 0;
  std::vector<std::int64_t> node_tags;
  std::vector<double> points;
  std::vector<std::int64_t> element_nodes;
  std::vector<std::int32_t> element_physical_tags;
  std::vector<double> element_metrics;
  std::vector<std::int64_t> face_nodes;
  std::vector<std::int32_t> face_elements;
  std::vector<std::int32_t> face_local_edges;
  std::vector<double> face_normals;
  std::vector<double> face_measures;
  std::vector<std::int32_t> face_physical_tags;
  std::vector<std::int32_t> physical_name_dimensions;
  std::vector<std::int32_t> physical_name_tags;
  std::vector<std::uint64_t> physical_name_offsets;
  std::vector<unsigned char> physical_name_bytes;
  // Partition metadata (ARCHITECTURE.md). The defaults describe an
  // unpartitioned mesh. The ghost layer is derived on demand from ownership and
  // topology, never persisted.
  int partition_rank_count = 1;
  std::vector<std::int32_t> element_owners;
  std::vector<std::int32_t> face_owners;

  std::size_t element_count() const noexcept { return element_physical_tags.size(); }
  std::size_t face_count() const noexcept { return face_physical_tags.size(); }

  // Line: Jacobian. Triangle: J00, J01, J10, J11, det(J).
  std::size_t metrics_per_element() const noexcept { return dimension == 1 ? 1U : 5U; }
};

} // namespace wasl
