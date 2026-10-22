#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "wasl/io/checkpoint.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/history_csv.hpp"
#include "wasl/io/prepared_mesh_file.hpp"
#include "wasl/io/pvd_writer.hpp"
#include "wasl/io/pvtu_writer.hpp"
#include "wasl/io/vtu_writer.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/partition/partition.hpp"
#include "wasl/mesh/prepared/owned_mesh_slice.hpp"
#include "wasl/mesh/prepared/prepare_mesh.hpp"
#include "wasl/mesh/prepared/prepared_mesh.hpp"

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string read_text(const std::filesystem::path &path) {
  std::ifstream input(path);
  std::ostringstream text;
  text << input.rdbuf();
  return text.str();
}

} // namespace

int main() {
  try {
    const auto directory = std::filesystem::temp_directory_path() / "wasl-case-io-test";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);

    wasl::LineMesh source;
    source.elements = {{10, 20, 0.0, 0.5, 7}, {20, 30, 0.5, 1.0, 7}};
    source.faces = {{0, -1, 0.0, -1.0, 3}, {0, 1, 0.5, 1.0, -1}, {1, -1, 1.0, 1.0, 4}};
    source.physical_names = {{{0, 3}, "left"}, {{0, 4}, "right"}};
    const auto prepared = wasl::prepare_mesh(source);
    require(prepared.element_nodes == std::vector<std::int64_t>({0, 1, 1, 2}) &&
                prepared.element_metrics == std::vector<double>({0.25, 0.25}) &&
                prepared.physical_name_offsets == std::vector<std::uint64_t>({0, 4, 9}),
            "line mesh preparation failed");

    wasl::PreparedMesh mesh{1,
                            2,
                            {10, 20, 30},
                            {0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0, 1, 1, 2},
                            {7, 7},
                            {0.25, 0.25},
                            {0, 0, 1, 1, 2, 2},
                            {0, -1, 0, 1, 1, -1},
                            {0, 0, 1, 0, 1, 0},
                            {-1.0, 1.0, 1.0},
                            {1.0, 1.0, 1.0},
                            {3, -1, 4},
                            {0, 0},
                            {3, 4},
                            {0, 4, 9},
                            {'l', 'e', 'f', 't', 'r', 'i', 'g', 'h', 't'},
                            1,
                            {0, 0},
                            {0, 0, 0}};
    wasl::write_prepared_mesh(directory / "prepared.h5", mesh);
    const auto loaded_mesh = wasl::read_prepared_mesh(directory / "prepared.h5");
    require(loaded_mesh.node_tags == mesh.node_tags &&
                loaded_mesh.face_elements == mesh.face_elements &&
                loaded_mesh.physical_name_bytes == mesh.physical_name_bytes,
            "prepared mesh round trip failed");
    require(loaded_mesh.partition_rank_count == 1 &&
                loaded_mesh.element_owners == std::vector<std::int32_t>({0, 0}) &&
                loaded_mesh.face_owners == std::vector<std::int32_t>({0, 0, 0}),
            "prepared mesh round trip lost partition ownership");

    // A version-1 prepared mesh carries no partition datasets. The reader must
    // back-fill a single-rank partition so callers see uniform metadata. Build a
    // genuine v1 file by writing v2, then downgrading the version attribute and
    // deleting the partition datasets.
    wasl::write_prepared_mesh(directory / "v1.h5", mesh);
    {
      H5::H5File file(directory / "v1.h5", H5F_ACC_RDWR);
      file.removeAttr("format_version");
      H5::DataSpace scalar(H5S_SCALAR);
      const int one = 1;
      file.createAttribute("format_version", H5::PredType::NATIVE_INT, scalar)
          .write(H5::PredType::NATIVE_INT, &one);
      file.unlink("element_owners");
      file.unlink("face_owners");
    }
    const auto legacy = wasl::read_prepared_mesh(directory / "v1.h5");
    require(legacy.partition_rank_count == 1 &&
                legacy.element_owners == std::vector<std::int32_t>({0, 0}) &&
                legacy.face_owners == std::vector<std::int32_t>({0, 0, 0}) &&
                legacy.face_elements == mesh.face_elements,
            "version-1 prepared mesh did not back-fill a single-rank partition");

    // A multi-rank partition must round-trip through HDF5 and pass validation.
    // This exercises the density check and the face-owner-equals-left rule at
    // rank_count > 1, which the single-rank path never reaches.
    wasl::PreparedMesh split = mesh;
    const auto partition = wasl::partition_mesh(split, 2);
    split.partition_rank_count = partition.rank_count;
    split.element_owners = partition.element_owners;
    split.face_owners = partition.face_owners;
    wasl::write_prepared_mesh(directory / "split.h5", split);
    const auto loaded_split = wasl::read_prepared_mesh(directory / "split.h5");
    require(loaded_split.partition_rank_count == 2 &&
                loaded_split.element_owners == std::vector<std::int32_t>({0, 1}) &&
                loaded_split.face_owners == std::vector<std::int32_t>({0, 0, 1}),
            "multi-rank partition did not round-trip");
    const auto rank1_mesh = wasl::owned_mesh_slice(loaded_split, 1);
    require(rank1_mesh.node_tags == std::vector<std::int64_t>({20, 30}) &&
                rank1_mesh.element_nodes == std::vector<std::int64_t>({0, 1}),
            "owned mesh slice did not compact local node indices");

    const wasl::Checkpoint checkpoint{0.25, 12, 1, 2, 2, "prepared.h5", {1.0, 2.0, 3.0, 4.0},
                                      0,    1,  0, {}};
    wasl::write_checkpoint(directory / "checkpoint.h5", checkpoint);
    const auto loaded_checkpoint = wasl::read_checkpoint(directory / "checkpoint.h5");
    require(loaded_checkpoint.time == checkpoint.time &&
                loaded_checkpoint.coefficients == checkpoint.coefficients &&
                loaded_checkpoint.prepared_mesh == checkpoint.prepared_mesh,
            "checkpoint round trip failed");
    wasl::write_checkpoint(directory / "checkpoint-v1.h5", checkpoint);
    {
      H5::H5File file(directory / "checkpoint-v1.h5", H5F_ACC_RDWR);
      file.removeAttr("format_version");
      H5::DataSpace scalar(H5S_SCALAR);
      const int one = 1;
      file.createAttribute("format_version", H5::PredType::NATIVE_INT, scalar)
          .write(H5::PredType::NATIVE_INT, &one);
      file.removeAttr("partition_rank");
      file.removeAttr("partition_rank_count");
      file.removeAttr("global_element_count");
      file.unlink("element_indices");
    }
    const auto legacy_checkpoint = wasl::read_checkpoint(directory / "checkpoint-v1.h5");
    require(legacy_checkpoint.coefficients == checkpoint.coefficients &&
                legacy_checkpoint.partition_rank_count == 1,
            "version-1 checkpoint did not back-fill serial ownership");

    wasl::write_csv_history(directory / "history.csv", {{0, 0.0, 0.1, 2.0}, {1, 0.1, 0.1, 2.0}});
    require(read_text(directory / "history.csv") ==
                "step,time,dt,mass\n0,0,0.10000000000000001,2\n1,0.10000000000000001,0."
                "10000000000000001,2\n",
            "CSV output is not deterministic");

    wasl::write_vtu(directory / "state.vtu", mesh, {1.0, 2.0});
    wasl::write_vtu(directory / "rank-1.vtu", rank1_mesh, {2.0});
    wasl::write_pvtu(directory / "state.pvtu", {"rank-0.vtu", "rank-1.vtu"});
    wasl::write_pvd(directory / "state.pvd", {{0.25, "state.vtu"}});
    require(read_text(directory / "state.vtu").find("NumberOfCells=\"2\"") != std::string::npos,
            "VTU output missing cell count");
    require(read_text(directory / "state.pvd").find("file=\"state.vtu\"") != std::string::npos,
            "PVD output missing snapshot");
    require(read_text(directory / "state.pvtu").find("rank-1.vtu") != std::string::npos,
            "PVTU output missing rank piece");

    std::filesystem::remove_all(directory);
    std::cout << "case I/O checks passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
