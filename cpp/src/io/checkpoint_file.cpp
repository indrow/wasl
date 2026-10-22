#include "wasl/io/checkpoint_file.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include "io/hdf5_io.hpp"
#include "io/text_output.hpp"

namespace wasl {
namespace {

constexpr int checkpoint_format_version = 2;

} // namespace

void write_checkpoint(const std::filesystem::path &path, const Checkpoint &checkpoint) {
  validate_checkpoint(checkpoint);
  create_parent_directories(path);
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
  if (!std::filesystem::is_regular_file(path))
    throw std::runtime_error("checkpoint HDF5 read: file does not exist: " + path.string());
  H5::Exception::dontPrint();
  try {
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    const int version = read_attribute<int>(file, "format_version");
    if (version < earliest_supported_format_version || version > checkpoint_format_version)
      throw std::runtime_error("checkpoint HDF5 read: unsupported format_version");
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
    validate_checkpoint(checkpoint);
    return checkpoint;
  } catch (const H5::Exception &error) {
    throw std::runtime_error("checkpoint HDF5 read: " + error.getDetailMsg());
  }
}

} // namespace wasl
