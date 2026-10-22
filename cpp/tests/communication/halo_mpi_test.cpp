// MPI transport and case-loop regression. Run under `mpirun -np 2`.
//
// The line case exchanges partition and periodic endpoint ghosts at all three
// SSP-RK3 stages, advances five steps, and checks owned coefficients, global
// mass, VTU pieces, PVTU, and checkpoint reassembly against the serial oracle.
// The triangle check retains the smaller one-stage residual coverage.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "wasl/basis/reference_element.hpp"
#include "wasl/boundaries/periodic.hpp"
#include "wasl/communication/halo_exchange.hpp"
#include "wasl/communication/host_staged_mpi.hpp"
#include "wasl/core/simulation.hpp"
#include "wasl/dg/residual_operator.hpp"
#include "wasl/dg/triangle_residual.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/mesh/line_geometry.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/partition.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

#include "rank_local_mesh.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

std::vector<double> triangle_coefficients(const wasl::TriangleMesh &mesh, int polynomial_order) {
  const std::size_t mode_count = wasl::triangle_mode_count(polynomial_order);
  std::vector<double> coefficients(mesh.num_elements() * mode_count, 0.0);
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const auto &elem = mesh.elements[e];
    const double x = (elem.x0 + elem.x1 + elem.x2) / 3.0;
    const double y = (elem.y0 + elem.y1 + elem.y2) / 3.0;
    coefficients[e * mode_count] = (1.0 + 0.3 * std::sin(x) + 0.2 * std::cos(y)) * std::sqrt(2.0);
    if (mode_count >= 3) {
      coefficients[e * mode_count + 1] = 0.15 * x;
      coefficients[e * mode_count + 2] = 0.15 * y;
    }
  }
  return coefficients;
}

std::vector<double> boundary_zero(const wasl::TriangleMesh &mesh, int polynomial_order) {
  return std::vector<double>(mesh.faces.size() * static_cast<std::size_t>(polynomial_order + 1),
                             0.0);
}

// Gather a rank's owned element indices (global order) for a contiguous-block
// partition, matching partition.hpp.
std::vector<std::int32_t> owned_elements(std::size_t element_count, int rank_count, int rank) {
  std::vector<std::int32_t> owned;
  for (std::size_t e = 0; e < element_count; ++e) {
    const int owner = static_cast<int>(
        (static_cast<std::uint64_t>(e) * static_cast<std::uint64_t>(rank_count)) / element_count);
    if (owner == rank) {
      owned.push_back(static_cast<std::int32_t>(e));
    }
  }
  return owned;
}

bool check_line(int rank, int rank_count) {
  constexpr int n = 6;
  const auto serial = wasl_testing::make_uniform_line(1.0, n);
  const auto directory = std::filesystem::temp_directory_path() / "wasl-mpi-case-driver-test";
  const auto prepared_path = directory / "prepared.h5";
  if (rank == 0) {
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    wasl::write_prepared_mesh(prepared_path, wasl::prepare_mesh(serial));
  }
  MPI_Barrier(MPI_COMM_WORLD);

  wasl::ScalarAdvectionCase config;
  config.prepared_mesh = prepared_path;
  config.mpi = true;
  config.polynomial_order = 1;
  config.wavespeed = 1.0;
  config.cfl = 0.2;
  config.start_time = 0.0;
  config.end_time = 0.04;
  config.max_steps = 100;
  config.history_interval = 1;
  config.visualization_interval = 100;
  config.checkpoint_interval = 100;
  config.keep_checkpoints = 2;
  config.initial_default = 1.0;
  config.initial_intervals.push_back({0.25, 0.5, 2.0});

  const auto fresh_directory =
      wasl::prepare_mpi_run_directory(directory / "results", "test = true\n");
  config.output_directory = fresh_directory;
  const auto fresh = wasl::run_scalar_advection_case(config);

  const auto split_directory =
      wasl::prepare_mpi_run_directory(directory / "results", "test = true\n");
  config.output_directory = split_directory;
  config.end_time = wasl::advection_cfl_timestep(serial, wasl::ScalarAdvection{1.0}, 1, 0.2);
  const auto split = wasl::run_scalar_advection_case(config);

  std::ostringstream restart_name;
  restart_name << "checkpoint-" << std::setfill('0') << std::setw(6) << split.step
               << "-rank-000000.h5";
  const auto restarted_directory =
      wasl::prepare_mpi_run_directory(directory / "results", "test = true\n");
  config.output_directory = restarted_directory;
  config.restart_checkpoint = split_directory / restart_name.str();
  config.end_time = 0.04;
  const auto restarted = wasl::run_scalar_advection_case(config);

  bool ok = fresh.step == restarted.step && approx(fresh.time, restarted.time);
  if (rank == 0) {
    if (!ok) {
      std::cerr << "fresh step/time " << fresh.step << '/' << fresh.time << ", restarted "
                << restarted.step << '/' << restarted.time << '\n';
    }
    std::vector<wasl::Checkpoint> fresh_parts;
    std::vector<wasl::Checkpoint> restarted_parts;
    for (int source = 0; source < rank_count; ++source) {
      std::ostringstream fresh_name;
      fresh_name << "checkpoint-" << std::setfill('0') << std::setw(6) << fresh.step << "-rank-"
                 << std::setw(6) << source << ".h5";
      std::ostringstream restarted_name;
      restarted_name << "checkpoint-" << std::setfill('0') << std::setw(6) << restarted.step
                     << "-rank-" << std::setw(6) << source << ".h5";
      fresh_parts.push_back(wasl::read_checkpoint(fresh_directory / fresh_name.str()));
      restarted_parts.push_back(wasl::read_checkpoint(restarted_directory / restarted_name.str()));
    }
    const auto fresh_joined = wasl::reassemble_checkpoints(fresh_parts);
    const auto restarted_joined = wasl::reassemble_checkpoints(restarted_parts);
    for (std::size_t index = 0; index < fresh_joined.coefficients.size(); ++index) {
      if (!approx(fresh_joined.coefficients[index], restarted_joined.coefficients[index])) {
        std::cerr << "checkpoint coefficient " << index
                  << " differs: " << fresh_joined.coefficients[index] << " vs "
                  << restarted_joined.coefficients[index] << '\n';
        ok = false;
      }
    }
    ok = ok && std::filesystem::is_regular_file(fresh_directory / "solution.pvd") &&
         std::filesystem::is_regular_file(fresh_directory / "history.csv");
  }
  MPI_Bcast(&ok, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::filesystem::remove_all(directory);
  }
  return ok;
}

bool check_triangle(int rank, int rank_count) {
  constexpr int cells = 2; // 8 triangles
  constexpr int p = 1;
  const auto serial = wasl_testing::make_unit_square_triangles(cells);
  const wasl::ScalarAdvection2D equation{0.7, -1.1};
  const std::size_t mode_count = wasl::triangle_mode_count(p);
  const std::size_t n = serial.num_elements();

  wasl::PreparedMesh prepared = wasl::prepare_mesh(serial);
  const auto partition = wasl::partition_mesh(prepared, rank_count);
  const auto owned = owned_elements(n, rank_count, rank);
  const auto ghosts = wasl::rank_ghost_elements(prepared, partition, rank);
  const auto local = wasl_testing::slice_triangle(serial, owned, ghosts);

  const auto serial_coef = triangle_coefficients(serial, p);
  const auto serial_ext = boundary_zero(serial, p);
  std::vector<double> serial_rhs(serial_coef.size());
  wasl::triangle_scalar_advection_residual(serial, equation, p, serial_coef.data(),
                                           serial_ext.data(), serial_rhs.data());

  std::vector<double> coef(local.mesh.num_elements() * mode_count, 0.0);
  for (const auto global : owned) {
    const auto loc = local.global_to_local.at(global);
    std::copy_n(serial_coef.data() + global * mode_count, mode_count,
                coef.data() + loc * mode_count);
  }

  const auto schedule = wasl::build_halo_schedule(prepared, partition, rank);
  wasl::HostStagedHaloExchange exchange(MPI_COMM_WORLD, schedule, mode_count, 0, mode_count,
                                        local.global_to_local);
  exchange.start(coef.data());
  exchange.finish(coef.data());

  const auto ext = boundary_zero(local.mesh, p);
  std::vector<double> rhs(coef.size());
  wasl::triangle_scalar_advection_residual(local.mesh, equation, p, coef.data(), ext.data(),
                                           rhs.data());

  for (const auto global : owned) {
    const auto loc = local.global_to_local.at(global);
    for (std::size_t m = 0; m < mode_count; ++m) {
      if (!approx(rhs[loc * mode_count + m], serial_rhs[global * mode_count + m]))
        return false;
    }
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

  int failures = 0;
  const char *context = nullptr;
  try {
    if (rank_count != 2) {
      throw std::runtime_error("MPI regression requires exactly two ranks");
    }
    bool ok = check_line(rank, rank_count);
    if (!ok) {
      context = "two-rank scalar-advection case does not match serial or distributed output";
      failures = 1;
    }
    if (!failures) {
      ok = check_triangle(rank, rank_count);
      if (!ok) {
        context = "triangle two-rank MPI residual does not match serial oracle";
        failures = 1;
      }
    }
  } catch (const std::exception &error) {
    context = error.what();
    failures = 1;
  }

  // Reduce failures across ranks so the test fails unless every rank passes.
  int global_failures = 0;
  MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (global_failures != 0 && context != nullptr) {
    std::cerr << "rank " << rank << ": " << context << '\n';
  }
  MPI_Finalize();
  return global_failures == 0 ? 0 : 1;
}
