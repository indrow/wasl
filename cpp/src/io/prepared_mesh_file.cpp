#include "wasl/io/prepared_mesh_file.hpp"

#include <stdexcept>
#include <string>

#include "io/hdf5_io.hpp"
#include "io/text_output.hpp"
#include "wasl/mesh/partition/partition.hpp"
#include "wasl/mesh/prepared/prepared_mesh_validation.hpp"

namespace wasl {
namespace {

constexpr int prepared_mesh_format_version = 2;

void read_partition(const H5::H5File &file, int version, PreparedMesh &mesh) {
  if (version >= 2) {
    mesh.partition_rank_count = read_attribute<int>(file, "partition_rank_count");
    mesh.element_owners = read_vector<std::int32_t>(file, "element_owners");
    mesh.face_owners = read_vector<std::int32_t>(file, "face_owners");
    return;
  }
  // Version 1 has no ownership; recompute a single-rank partition, checking
  // the face array shape partition_mesh indexes first.
  if (mesh.face_elements.size() != mesh.face_count() * 2)
    throw std::runtime_error("prepared mesh: face array size mismatch");
  assign_partition(mesh, partition_mesh(mesh, 1));
}

} // namespace

void write_prepared_mesh(const std::filesystem::path &path, const PreparedMesh &mesh) {
  validate_prepared_mesh(mesh);
  create_parent_directories(path);
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_TRUNC);
    write_attribute(file, "format_version", prepared_mesh_format_version);
    write_attribute(file, "dimension", mesh.dimension);
    write_attribute(file, "nodes_per_element", mesh.nodes_per_element);
    write_vector(file, "node_tags", mesh.node_tags);
    write_vector(file, "points", mesh.points);
    write_vector(file, "element_nodes", mesh.element_nodes);
    write_vector(file, "element_physical_tags", mesh.element_physical_tags);
    write_vector(file, "element_metrics", mesh.element_metrics);
    write_vector(file, "face_nodes", mesh.face_nodes);
    write_vector(file, "face_elements", mesh.face_elements);
    write_vector(file, "face_local_edges", mesh.face_local_edges);
    write_vector(file, "face_normals", mesh.face_normals);
    write_vector(file, "face_measures", mesh.face_measures);
    write_vector(file, "face_physical_tags", mesh.face_physical_tags);
    write_vector(file, "physical_name_dimensions", mesh.physical_name_dimensions);
    write_vector(file, "physical_name_tags", mesh.physical_name_tags);
    write_vector(file, "physical_name_offsets", mesh.physical_name_offsets);
    write_vector(file, "physical_name_bytes", mesh.physical_name_bytes);
    write_attribute(file, "partition_rank_count", mesh.partition_rank_count);
    write_vector(file, "element_owners", mesh.element_owners);
    write_vector(file, "face_owners", mesh.face_owners);
  } catch (const H5::Exception &error) {
    throw std::runtime_error("prepared mesh HDF5 write: " + error.getDetailMsg());
  }
}

PreparedMesh read_prepared_mesh(const std::filesystem::path &path) {
  if (!std::filesystem::is_regular_file(path))
    throw std::runtime_error("prepared mesh HDF5 read: file does not exist: " + path.string());
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    const int version = read_attribute<int>(file, "format_version");
    if (version < earliest_supported_format_version || version > prepared_mesh_format_version)
      throw std::runtime_error("prepared mesh HDF5 read: unsupported format_version");
    PreparedMesh mesh;
    mesh.dimension = read_attribute<int>(file, "dimension");
    mesh.nodes_per_element = read_attribute<int>(file, "nodes_per_element");
    mesh.node_tags = read_vector<std::int64_t>(file, "node_tags");
    mesh.points = read_vector<double>(file, "points");
    mesh.element_nodes = read_vector<std::int64_t>(file, "element_nodes");
    mesh.element_physical_tags = read_vector<std::int32_t>(file, "element_physical_tags");
    mesh.element_metrics = read_vector<double>(file, "element_metrics");
    mesh.face_nodes = read_vector<std::int64_t>(file, "face_nodes");
    mesh.face_elements = read_vector<std::int32_t>(file, "face_elements");
    mesh.face_local_edges = read_vector<std::int32_t>(file, "face_local_edges");
    mesh.face_normals = read_vector<double>(file, "face_normals");
    mesh.face_measures = read_vector<double>(file, "face_measures");
    mesh.face_physical_tags = read_vector<std::int32_t>(file, "face_physical_tags");
    mesh.physical_name_dimensions = read_vector<std::int32_t>(file, "physical_name_dimensions");
    mesh.physical_name_tags = read_vector<std::int32_t>(file, "physical_name_tags");
    mesh.physical_name_offsets = read_vector<std::uint64_t>(file, "physical_name_offsets");
    mesh.physical_name_bytes = read_vector<unsigned char>(file, "physical_name_bytes");
    read_partition(file, version, mesh);
    validate_prepared_mesh(mesh);
    return mesh;
  } catch (const H5::Exception &error) {
    throw std::runtime_error("prepared mesh HDF5 read: " + error.getDetailMsg());
  }
}

} // namespace wasl
