#include "wasl/io/case_io.hpp"

#include <H5Cpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/partition.hpp"
#include "wasl/mesh/triangle_mesh.hpp"

namespace wasl {
namespace {

// Prepared-mesh format history:
//   1: connectivity, metrics, faces, physical names.
//   2: adds partition ownership (rank_count, element_owners, face_owners).
// Reads of version 1 back-fill a single-rank partition; writes always emit 2.
constexpr int prepared_mesh_format_version = 2;
constexpr int checkpoint_format_version = 2;
constexpr int earliest_supported_format_version = 1;

template <typename T> const H5::PredType &h5_type();
template <> const H5::PredType &h5_type<double>() { return H5::PredType::NATIVE_DOUBLE; }
template <> const H5::PredType &h5_type<std::int64_t>() { return H5::PredType::NATIVE_INT64; }
template <> const H5::PredType &h5_type<std::int32_t>() { return H5::PredType::NATIVE_INT32; }
template <> const H5::PredType &h5_type<std::uint64_t>() { return H5::PredType::NATIVE_UINT64; }
template <> const H5::PredType &h5_type<unsigned char>() { return H5::PredType::NATIVE_UCHAR; }

template <typename T> void write_attribute(H5::H5Object &object, const char *name, const T &value) {
  H5::DataSpace scalar(H5S_SCALAR);
  object.createAttribute(name, h5_type<T>(), scalar).write(h5_type<T>(), &value);
}

template <typename T> T read_attribute(const H5::H5Object &object, const char *name) {
  T value{};
  object.openAttribute(name).read(h5_type<T>(), &value);
  return value;
}

template <typename T>
void write_vector(H5::H5File &file, const char *name, const std::vector<T> &values) {
  const hsize_t size = values.size();
  H5::DataSpace space(1, &size);
  auto dataset = file.createDataSet(name, h5_type<T>(), space);
  if (!values.empty()) {
    dataset.write(values.data(), h5_type<T>());
  }
}

template <typename T> std::vector<T> read_vector(const H5::H5File &file, const char *name) {
  const auto dataset = file.openDataSet(name);
  const auto space = dataset.getSpace();
  if (space.getSimpleExtentNdims() != 1) {
    throw std::runtime_error(std::string("HDF5 I/O: dataset ") + name + " must be one-dimensional");
  }
  hsize_t size = 0;
  space.getSimpleExtentDims(&size);
  std::vector<T> values(static_cast<std::size_t>(size));
  if (!values.empty()) {
    dataset.read(values.data(), h5_type<T>());
  }
  return values;
}

void validate_prepared_mesh(const PreparedMesh &mesh) {
  if (mesh.dimension < 1 || mesh.dimension > 2) {
    throw std::runtime_error("prepared mesh: dimension must be 1 or 2");
  }
  if (mesh.nodes_per_element != mesh.dimension + 1) {
    throw std::runtime_error("prepared mesh: nodes_per_element does not match dimension");
  }
  if (mesh.points.size() != mesh.node_tags.size() * 3) {
    throw std::runtime_error("prepared mesh: points size does not match node tags");
  }
  if (std::any_of(mesh.points.begin(), mesh.points.end(),
                  [](const auto value) { return !std::isfinite(value); })) {
    throw std::runtime_error("prepared mesh: point coordinates must be finite");
  }
  if (mesh.element_nodes.empty() ||
      mesh.element_nodes.size() % static_cast<std::size_t>(mesh.nodes_per_element) != 0) {
    throw std::runtime_error("prepared mesh: invalid element connectivity");
  }
  const auto element_count = mesh.element_nodes.size() / mesh.nodes_per_element;
  if (std::any_of(mesh.element_nodes.begin(), mesh.element_nodes.end(), [&](const auto node) {
        return node < 0 || static_cast<std::size_t>(node) >= mesh.node_tags.size();
      })) {
    throw std::runtime_error("prepared mesh: element has invalid local node index");
  }
  if (mesh.element_physical_tags.size() != element_count) {
    throw std::runtime_error("prepared mesh: element physical-tag count mismatch");
  }
  const auto metrics_per_element = mesh.dimension == 1 ? 1U : 5U;
  if (mesh.element_metrics.size() != element_count * metrics_per_element) {
    throw std::runtime_error("prepared mesh: element metric count mismatch");
  }
  if (std::any_of(mesh.element_metrics.begin(), mesh.element_metrics.end(),
                  [](const auto value) { return !std::isfinite(value); })) {
    throw std::runtime_error("prepared mesh: element metrics must be finite");
  }
  const auto face_count = mesh.face_physical_tags.size();
  if (mesh.face_nodes.size() != face_count * 2 || mesh.face_elements.size() != face_count * 2 ||
      mesh.face_local_edges.size() != face_count * 2 ||
      mesh.face_normals.size() != face_count * static_cast<std::size_t>(mesh.dimension) ||
      mesh.face_measures.size() != face_count) {
    throw std::runtime_error("prepared mesh: face array size mismatch");
  }
  for (std::size_t face = 0; face < face_count; ++face) {
    const auto left = mesh.face_elements[face * 2];
    const auto right = mesh.face_elements[face * 2 + 1];
    if (left < 0 || static_cast<std::size_t>(left) >= element_count ||
        (right != invalid_element &&
         (right < 0 || static_cast<std::size_t>(right) >= element_count))) {
      throw std::runtime_error("prepared mesh: face has invalid element index");
    }
  }
  if (std::any_of(mesh.face_nodes.begin(), mesh.face_nodes.end(), [&](const auto node) {
        return node < 0 || static_cast<std::size_t>(node) >= mesh.node_tags.size();
      })) {
    throw std::runtime_error("prepared mesh: face has invalid local node index");
  }
  if (std::any_of(mesh.face_normals.begin(), mesh.face_normals.end(),
                  [](const auto value) { return !std::isfinite(value); }) ||
      std::any_of(mesh.face_measures.begin(), mesh.face_measures.end(),
                  [](const auto value) { return !std::isfinite(value) || value <= 0.0; })) {
    throw std::runtime_error("prepared mesh: face geometry must be finite with positive measures");
  }
  if (mesh.physical_name_offsets.empty() ||
      mesh.physical_name_dimensions.size() != mesh.physical_name_tags.size() ||
      mesh.physical_name_offsets.size() != mesh.physical_name_tags.size() + 1 ||
      mesh.physical_name_offsets.front() != 0 ||
      mesh.physical_name_offsets.back() != mesh.physical_name_bytes.size() ||
      !std::is_sorted(mesh.physical_name_offsets.begin(), mesh.physical_name_offsets.end())) {
    throw std::runtime_error("prepared mesh: physical-name array size mismatch");
  }
  if (mesh.partition_rank_count < 1) {
    throw std::runtime_error("prepared mesh: partition_rank_count must be positive");
  }
  if (mesh.element_owners.size() != element_count) {
    throw std::runtime_error("prepared mesh: element owner count mismatch");
  }
  if (mesh.face_owners.size() != face_count) {
    throw std::runtime_error("prepared mesh: face owner count mismatch");
  }
  // Every owner must name a real rank.
  for (const auto owner : mesh.element_owners) {
    if (owner < 0 || owner >= mesh.partition_rank_count) {
      throw std::runtime_error("prepared mesh: element owner is out of range");
    }
  }
  for (const auto owner : mesh.face_owners) {
    if (owner < 0 || owner >= mesh.partition_rank_count) {
      throw std::runtime_error("prepared mesh: face owner is out of range");
    }
  }
  // When there are at least as many elements as ranks, the contiguous block
  // scheme assigns at least one element to every rank, so each rank gets work.
  // With more ranks than elements some ranks are legitimately idle.
  if (element_count >= static_cast<std::size_t>(mesh.partition_rank_count)) {
    std::vector<unsigned char> owns_rank(mesh.partition_rank_count, 0);
    for (const auto owner : mesh.element_owners) {
      owns_rank[static_cast<std::size_t>(owner)] = 1;
    }
    for (std::size_t r = 0; r < static_cast<std::size_t>(mesh.partition_rank_count); ++r) {
      if (!owns_rank[r]) {
        throw std::runtime_error("prepared mesh: a rank owns no elements");
      }
    }
  }
  // A face is owned by its left element's owner (NUMERICAL_CONVENTIONS.md).
  for (std::size_t face = 0; face < face_count; ++face) {
    const auto left = mesh.face_elements[face * 2];
    const auto expected = mesh.element_owners.at(static_cast<std::size_t>(left));
    if (mesh.face_owners[face] != expected) {
      throw std::runtime_error("prepared mesh: face owner does not match its left element owner");
    }
  }
}

std::int32_t add_node(PreparedMesh &result, std::map<std::int64_t, std::int32_t> &indices,
                      std::int64_t tag, double x, double y) {
  const auto found = indices.find(tag);
  if (found != indices.end()) {
    return found->second;
  }
  if (result.node_tags.size() >
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::runtime_error("prepared mesh: local node count exceeds int32 range");
  }
  const auto index = static_cast<std::int32_t>(result.node_tags.size());
  indices.emplace(tag, index);
  result.node_tags.push_back(tag);
  result.points.insert(result.points.end(), {x, y, 0.0});
  return index;
}

void require_parent(const std::filesystem::path &path) {
  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }
}

std::ofstream open_text(const std::filesystem::path &path) {
  require_parent(path);
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("file output: cannot open " + path.string());
  }
  output << std::setprecision(17);
  return output;
}

template <typename Mesh> void add_physical_names(PreparedMesh &result, const Mesh &mesh) {
  result.physical_name_offsets.push_back(0);
  for (const auto &[key, name] : mesh.physical_names) {
    result.physical_name_dimensions.push_back(key.first);
    result.physical_name_tags.push_back(key.second);
    result.physical_name_bytes.insert(result.physical_name_bytes.end(), name.begin(), name.end());
    result.physical_name_offsets.push_back(result.physical_name_bytes.size());
  }
}

std::string xml_escape(const std::string &value) {
  std::string result;
  for (const char character : value) {
    switch (character) {
    case '&':
      result += "&amp;";
      break;
    case '<':
      result += "&lt;";
      break;
    case '"':
      result += "&quot;";
      break;
    default:
      result += character;
    }
  }
  return result;
}

bool checkpoint_layout_matches(const Checkpoint &checkpoint) {
  if (checkpoint.variable_count == 0 || checkpoint.mode_count == 0 ||
      checkpoint.element_count == 0 || checkpoint.prepared_mesh.empty() ||
      checkpoint.variable_count >
          std::numeric_limits<std::uint64_t>::max() / checkpoint.mode_count) {
    return false;
  }
  const auto variable_modes = checkpoint.variable_count * checkpoint.mode_count;
  return variable_modes <= std::numeric_limits<std::uint64_t>::max() / checkpoint.element_count &&
         variable_modes * checkpoint.element_count == checkpoint.coefficients.size();
}

std::uint64_t checkpoint_global_element_count(const Checkpoint &checkpoint) {
  return checkpoint.global_element_count == 0 ? checkpoint.element_count
                                              : checkpoint.global_element_count;
}

std::vector<std::int64_t> checkpoint_element_indices(const Checkpoint &checkpoint) {
  if (!checkpoint.element_indices.empty()) {
    return checkpoint.element_indices;
  }
  std::vector<std::int64_t> indices(checkpoint.element_count);
  std::iota(indices.begin(), indices.end(), 0);
  return indices;
}

void validate_checkpoint_partition(const Checkpoint &checkpoint) {
  if (!checkpoint_layout_matches(checkpoint)) {
    throw std::runtime_error("checkpoint: coefficient size does not match declared layout");
  }
  if (checkpoint.partition_rank_count < 1 || checkpoint.partition_rank < 0 ||
      checkpoint.partition_rank >= checkpoint.partition_rank_count) {
    throw std::runtime_error("checkpoint: invalid partition rank");
  }
  const auto global_count = checkpoint_global_element_count(checkpoint);
  const auto indices = checkpoint_element_indices(checkpoint);
  if (indices.size() != checkpoint.element_count || global_count < checkpoint.element_count) {
    throw std::runtime_error("checkpoint: invalid partition element count");
  }
  std::vector<unsigned char> seen(global_count, 0);
  for (const auto index : indices) {
    if (index < 0 || static_cast<std::uint64_t>(index) >= global_count ||
        seen[static_cast<std::size_t>(index)] != 0) {
      throw std::runtime_error("checkpoint: invalid partition element index");
    }
    seen[static_cast<std::size_t>(index)] = 1;
  }
}

} // namespace

PreparedMesh prepare_mesh(const LineMesh &mesh) {
  PreparedMesh result;
  result.dimension = 1;
  result.nodes_per_element = 2;
  std::map<std::int64_t, std::int32_t> indices;
  for (const auto &element : mesh.elements) {
    const auto left = add_node(result, indices, element.left_node, element.x_left, 0.0);
    const auto right = add_node(result, indices, element.right_node, element.x_right, 0.0);
    result.element_nodes.insert(result.element_nodes.end(), {left, right});
    result.element_physical_tags.push_back(element.physical_tag);
    result.element_metrics.push_back((element.x_right - element.x_left) / 2.0);
  }
  for (const auto &face : mesh.faces) {
    const auto &element = mesh.elements.at(static_cast<std::size_t>(face.left_element));
    const auto node = face.x == element.x_left ? element.left_node : element.right_node;
    result.face_nodes.insert(result.face_nodes.end(), {indices.at(node), indices.at(node)});
    result.face_elements.insert(result.face_elements.end(),
                                {face.left_element, face.right_element});
    result.face_local_edges.insert(result.face_local_edges.end(),
                                   {face.outward_normal < 0.0 ? 0 : 1, 0});
    result.face_normals.push_back(face.outward_normal);
    result.face_measures.push_back(1.0);
    result.face_physical_tags.push_back(face.physical_tag);
  }
  add_physical_names(result, mesh);
  // Default single-rank partition; multi-rank partitioning is applied later
  // from the bridge once a rank count is known (ARCHITECTURE.md).
  const auto partition = partition_mesh(result, 1);
  result.partition_rank_count = partition.rank_count;
  result.element_owners = std::move(partition.element_owners);
  result.face_owners = std::move(partition.face_owners);
  return result;
}

PreparedMesh prepare_mesh(const TriangleMesh &mesh) {
  PreparedMesh result;
  result.dimension = 2;
  result.nodes_per_element = 3;
  std::map<std::int64_t, std::int32_t> indices;
  for (const auto &element : mesh.elements) {
    const auto node0 = add_node(result, indices, element.node0, element.x0, element.y0);
    const auto node1 = add_node(result, indices, element.node1, element.x1, element.y1);
    const auto node2 = add_node(result, indices, element.node2, element.x2, element.y2);
    result.element_nodes.insert(result.element_nodes.end(), {node0, node1, node2});
    result.element_physical_tags.push_back(element.physical_tag);
    const auto map =
        TriangleMap{element.x0, element.y0, element.x1, element.y1, element.x2, element.y2};
    result.element_metrics.insert(result.element_metrics.end(),
                                  {map.ax() / 2.0, map.bx() / 2.0, map.ay() / 2.0, map.by() / 2.0,
                                   map.jacobian_determinant()});
  }
  for (const auto &face : mesh.faces) {
    result.face_nodes.insert(result.face_nodes.end(),
                             {indices.at(face.vertex_a), indices.at(face.vertex_b)});
    result.face_elements.insert(result.face_elements.end(),
                                {face.left_element, face.right_element});
    result.face_local_edges.insert(result.face_local_edges.end(),
                                   {face.left_local_edge, face.right_local_edge});
    result.face_normals.insert(result.face_normals.end(),
                               {face.outward_normal_x, face.outward_normal_y});
    const auto &owner = mesh.elements.at(static_cast<std::size_t>(face.left_element));
    const auto map = TriangleMap{owner.x0, owner.y0, owner.x1, owner.y1, owner.x2, owner.y2};
    result.face_measures.push_back(map.edge_length(face.left_local_edge));
    result.face_physical_tags.push_back(face.physical_tag);
  }
  add_physical_names(result, mesh);
  // Default single-rank partition; multi-rank partitioning is applied later
  // from the bridge once a rank count is known (ARCHITECTURE.md).
  const auto partition = partition_mesh(result, 1);
  result.partition_rank_count = partition.rank_count;
  result.element_owners = std::move(partition.element_owners);
  result.face_owners = std::move(partition.face_owners);
  return result;
}

void write_prepared_mesh(const std::filesystem::path &path, const PreparedMesh &mesh) {
  validate_prepared_mesh(mesh);
  require_parent(path);
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
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("prepared mesh HDF5 read: file does not exist: " + path.string());
  }
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    const int version = read_attribute<int>(file, "format_version");
    if (version < earliest_supported_format_version || version > prepared_mesh_format_version) {
      throw std::runtime_error("prepared mesh HDF5 read: unsupported format_version");
    }
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
    if (version >= 2) {
      mesh.partition_rank_count = read_attribute<int>(file, "partition_rank_count");
      mesh.element_owners = read_vector<std::int32_t>(file, "element_owners");
      mesh.face_owners = read_vector<std::int32_t>(file, "face_owners");
    } else {
      // Back-fill a single-rank partition for version-1 meshes so callers see
      // uniform metadata. Ownership is recomputed rather than read, so check
      // the face array shape partition_mesh indexes before calling it.
      if (mesh.face_elements.size() != mesh.face_physical_tags.size() * 2) {
        throw std::runtime_error("prepared mesh: face array size mismatch");
      }
      const auto partition = partition_mesh(mesh, 1);
      mesh.partition_rank_count = partition.rank_count;
      mesh.element_owners = partition.element_owners;
      mesh.face_owners = partition.face_owners;
    }
    validate_prepared_mesh(mesh);
    return mesh;
  } catch (const H5::Exception &error) {
    throw std::runtime_error("prepared mesh HDF5 read: " + error.getDetailMsg());
  }
}

void write_checkpoint(const std::filesystem::path &path, const Checkpoint &checkpoint) {
  validate_checkpoint_partition(checkpoint);
  require_parent(path);
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_TRUNC);
    write_attribute(file, "format_version", checkpoint_format_version);
    write_attribute(file, "time", checkpoint.time);
    write_attribute(file, "step", checkpoint.step);
    write_attribute(file, "variable_count", checkpoint.variable_count);
    write_attribute(file, "mode_count", checkpoint.mode_count);
    write_attribute(file, "element_count", checkpoint.element_count);
    write_vector(file, "prepared_mesh",
                 std::vector<unsigned char>(checkpoint.prepared_mesh.begin(),
                                            checkpoint.prepared_mesh.end()));
    write_vector(file, "coefficients", checkpoint.coefficients);
    write_attribute(file, "partition_rank", checkpoint.partition_rank);
    write_attribute(file, "partition_rank_count", checkpoint.partition_rank_count);
    write_attribute(file, "global_element_count", checkpoint_global_element_count(checkpoint));
    write_vector(file, "element_indices", checkpoint_element_indices(checkpoint));
  } catch (const H5::Exception &error) {
    throw std::runtime_error("checkpoint HDF5 write: " + error.getDetailMsg());
  }
}

Checkpoint read_checkpoint(const std::filesystem::path &path) {
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("checkpoint HDF5 read: file does not exist: " + path.string());
  }
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    const int version = read_attribute<int>(file, "format_version");
    if (version < earliest_supported_format_version || version > checkpoint_format_version) {
      throw std::runtime_error("checkpoint HDF5 read: unsupported format_version");
    }
    Checkpoint checkpoint;
    checkpoint.time = read_attribute<double>(file, "time");
    checkpoint.step = read_attribute<std::uint64_t>(file, "step");
    checkpoint.variable_count = read_attribute<std::uint64_t>(file, "variable_count");
    checkpoint.mode_count = read_attribute<std::uint64_t>(file, "mode_count");
    checkpoint.element_count = read_attribute<std::uint64_t>(file, "element_count");
    const auto mesh = read_vector<unsigned char>(file, "prepared_mesh");
    checkpoint.prepared_mesh.assign(mesh.begin(), mesh.end());
    checkpoint.coefficients = read_vector<double>(file, "coefficients");
    if (version >= 2) {
      checkpoint.partition_rank = read_attribute<int>(file, "partition_rank");
      checkpoint.partition_rank_count = read_attribute<int>(file, "partition_rank_count");
      checkpoint.global_element_count = read_attribute<std::uint64_t>(file, "global_element_count");
      checkpoint.element_indices = read_vector<std::int64_t>(file, "element_indices");
    }
    validate_checkpoint_partition(checkpoint);
    return checkpoint;
  } catch (const H5::Exception &error) {
    throw std::runtime_error("checkpoint HDF5 read: " + error.getDetailMsg());
  }
}

Checkpoint reassemble_checkpoints(const std::vector<Checkpoint> &parts) {
  if (parts.empty()) {
    throw std::runtime_error("checkpoint reassembly: no rank-local checkpoints");
  }
  const auto &first = parts.front();
  validate_checkpoint_partition(first);
  const auto global_count = checkpoint_global_element_count(first);
  if (parts.size() != static_cast<std::size_t>(first.partition_rank_count)) {
    throw std::runtime_error("checkpoint reassembly: rank count mismatch");
  }

  Checkpoint result;
  result.time = first.time;
  result.step = first.step;
  result.variable_count = first.variable_count;
  result.mode_count = first.mode_count;
  result.element_count = global_count;
  result.prepared_mesh = first.prepared_mesh;
  result.coefficients.resize(first.variable_count * first.mode_count * global_count);
  std::vector<unsigned char> seen_rank(first.partition_rank_count, 0);
  std::vector<unsigned char> seen_element(global_count, 0);
  for (const auto &part : parts) {
    validate_checkpoint_partition(part);
    if (part.time != first.time || part.step != first.step ||
        part.variable_count != first.variable_count || part.mode_count != first.mode_count ||
        checkpoint_global_element_count(part) != global_count ||
        part.partition_rank_count != first.partition_rank_count ||
        part.prepared_mesh != first.prepared_mesh ||
        seen_rank[static_cast<std::size_t>(part.partition_rank)] != 0) {
      throw std::runtime_error("checkpoint reassembly: inconsistent rank-local checkpoints");
    }
    seen_rank[static_cast<std::size_t>(part.partition_rank)] = 1;
    const auto indices = checkpoint_element_indices(part);
    for (std::size_t local = 0; local < indices.size(); ++local) {
      const auto global = static_cast<std::size_t>(indices[local]);
      if (seen_element[global] != 0) {
        throw std::runtime_error("checkpoint reassembly: duplicate owned element");
      }
      seen_element[global] = 1;
      for (std::uint64_t variable = 0; variable < part.variable_count; ++variable) {
        for (std::uint64_t mode = 0; mode < part.mode_count; ++mode) {
          const auto source = (variable * part.mode_count + mode) * part.element_count + local;
          const auto target = (variable * part.mode_count + mode) * global_count + global;
          result.coefficients[target] = part.coefficients[source];
        }
      }
    }
  }
  if (std::find(seen_element.begin(), seen_element.end(), 0) != seen_element.end()) {
    throw std::runtime_error("checkpoint reassembly: missing owned element");
  }
  return result;
}

void write_csv_history(const std::filesystem::path &path, const std::vector<HistoryRow> &rows) {
  auto output = open_text(path);
  output << "step,time,dt,mass\n";
  for (const auto &row : rows) {
    output << row.step << ',' << row.time << ',' << row.time_step << ',' << row.mass << '\n';
  }
}

void write_vtu(const std::filesystem::path &path, const PreparedMesh &mesh,
               const std::vector<double> &cell_values, std::size_t component_count) {
  validate_prepared_mesh(mesh);
  const auto cell_count = mesh.element_physical_tags.size();
  if (component_count == 0 || cell_values.size() != cell_count * component_count) {
    throw std::runtime_error("VTU output: cell-value size mismatch");
  }
  auto output = open_text(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<UnstructuredGrid><Piece NumberOfPoints=\""
         << mesh.node_tags.size() << "\" NumberOfCells=\"" << cell_count << "\">\n";
  output << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (double value : mesh.points) {
    output << value << ' ';
  }
  output << "\n</DataArray></Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" "
            "format=\"ascii\">\n";
  for (const auto index : mesh.element_nodes) {
    if (index < 0 || static_cast<std::size_t>(index) >= mesh.node_tags.size()) {
      throw std::runtime_error("VTU output: element references invalid local node index");
    }
    output << index << ' ';
  }
  output << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 1; i <= cell_count; ++i) {
    output << i * static_cast<std::size_t>(mesh.nodes_per_element) << ' ';
  }
  output << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  const int cell_type = mesh.dimension == 1 ? 3 : 5;
  for (std::size_t i = 0; i < cell_count; ++i) {
    output << cell_type << ' ';
  }
  output << "\n</DataArray>\n</Cells>\n<CellData Scalars=\"state\"><DataArray type=\"Float64\" "
            "Name=\"state\" NumberOfComponents=\""
         << component_count << "\" format=\"ascii\">\n";
  for (double value : cell_values) {
    output << value << ' ';
  }
  output << "\n</DataArray></CellData>\n</Piece></UnstructuredGrid></VTKFile>\n";
}

PreparedMesh owned_mesh_slice(const PreparedMesh &mesh, int rank) {
  validate_prepared_mesh(mesh);
  if (rank < 0 || rank >= mesh.partition_rank_count) {
    throw std::runtime_error("owned mesh slice: rank is out of range");
  }
  PreparedMesh result;
  result.dimension = mesh.dimension;
  result.nodes_per_element = mesh.nodes_per_element;
  result.physical_name_dimensions = mesh.physical_name_dimensions;
  result.physical_name_tags = mesh.physical_name_tags;
  result.physical_name_offsets = mesh.physical_name_offsets;
  result.physical_name_bytes = mesh.physical_name_bytes;
  std::map<std::int64_t, std::int32_t> indices;
  const auto metrics_per_element = mesh.dimension == 1 ? 1U : 5U;
  for (std::size_t element = 0; element < mesh.element_owners.size(); ++element) {
    if (mesh.element_owners[element] != rank) {
      continue;
    }
    for (int local_node = 0; local_node < mesh.nodes_per_element; ++local_node) {
      const auto source = static_cast<std::size_t>(
          mesh.element_nodes[element * static_cast<std::size_t>(mesh.nodes_per_element) +
                             static_cast<std::size_t>(local_node)]);
      result.element_nodes.push_back(add_node(result, indices, mesh.node_tags[source],
                                              mesh.points[source * 3],
                                              mesh.points[source * 3 + 1]));
    }
    result.element_physical_tags.push_back(mesh.element_physical_tags[element]);
    result.element_metrics.insert(
        result.element_metrics.end(), mesh.element_metrics.begin() + element * metrics_per_element,
        mesh.element_metrics.begin() + (element + 1) * metrics_per_element);
    result.element_owners.push_back(0);
  }
  if (result.element_nodes.empty()) {
    throw std::runtime_error("owned mesh slice: rank owns no elements");
  }
  return result;
}

void write_pvtu(const std::filesystem::path &path, const std::vector<std::string> &pieces,
                std::size_t component_count) {
  if (pieces.empty() || component_count == 0) {
    throw std::runtime_error("PVTU output: pieces and components must be non-empty");
  }
  auto output = open_text(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n<PUnstructuredGrid GhostLevel=\"0\">\n"
            "<PCellData Scalars=\"state\"><PDataArray type=\"Float64\" Name=\"state\" "
            "NumberOfComponents=\""
         << component_count
         << "\"/></PCellData>\n"
            "<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n";
  for (const auto &piece : pieces) {
    output << "<Piece Source=\"" << xml_escape(piece) << "\"/>\n";
  }
  output << "</PUnstructuredGrid></VTKFile>\n";
}

void write_pvd(const std::filesystem::path &path,
               const std::vector<std::pair<double, std::string>> &snapshots) {
  auto output = open_text(path);
  output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n<Collection>\n";
  for (const auto &[time, file] : snapshots) {
    output << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\""
           << xml_escape(file) << "\"/>\n";
  }
  output << "</Collection></VTKFile>\n";
}

} // namespace wasl
