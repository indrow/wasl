#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace wasl {

struct LineMesh;
struct TriangleMesh;

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
  // Partition metadata (ARCHITECTURE.md). Defaults describe an unpartitioned mesh,
  // so freshly constructed PreparedMesh values remain valid. The ghost layer
  // is derived on demand from ownership and topology, never persisted per rank.
  int partition_rank_count = 1;
  std::vector<std::int32_t> element_owners;
  std::vector<std::int32_t> face_owners;
};

PreparedMesh prepare_mesh(const LineMesh &mesh);
PreparedMesh prepare_mesh(const TriangleMesh &mesh);
void write_prepared_mesh(const std::filesystem::path &path, const PreparedMesh &mesh);
PreparedMesh read_prepared_mesh(const std::filesystem::path &path);

struct Checkpoint {
  double time = 0.0;
  std::uint64_t step = 0;
  std::uint64_t variable_count = 0;
  std::uint64_t mode_count = 0;
  std::uint64_t element_count = 0;
  std::string prepared_mesh;
  std::vector<double> coefficients;
  int partition_rank = 0;
  int partition_rank_count = 1;
  std::uint64_t global_element_count = 0;
  std::vector<std::int64_t> element_indices;
};

void write_checkpoint(const std::filesystem::path &path, const Checkpoint &checkpoint);
Checkpoint read_checkpoint(const std::filesystem::path &path);
Checkpoint reassemble_checkpoints(const std::vector<Checkpoint> &parts);

struct HistoryRow {
  std::uint64_t step;
  double time;
  double time_step;
  double mass;
};

void write_csv_history(const std::filesystem::path &path, const std::vector<HistoryRow> &rows);
void write_vtu(const std::filesystem::path &path, const PreparedMesh &mesh,
               const std::vector<double> &cell_values, std::size_t component_count = 1);
PreparedMesh owned_mesh_slice(const PreparedMesh &mesh, int rank);
void write_pvtu(const std::filesystem::path &path, const std::vector<std::string> &pieces,
                std::size_t component_count = 1);
void write_pvd(const std::filesystem::path &path,
               const std::vector<std::pair<double, std::string>> &snapshots);

} // namespace wasl
