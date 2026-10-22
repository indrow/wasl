#include "case/scalar_advection_drivers.hpp"

#include <stdexcept>

#ifdef WASL_ENABLE_MPI

#include <algorithm>
#include <cmath>
#include <deque>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include "case/line_boundary_resolution.hpp"
#include "communication/mpi_session.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/line_exterior.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/case/scalar_coefficient_layout.hpp"
#include "wasl/communication/halo_schedule.hpp"
#include "wasl/communication/host_staged_mpi.hpp"
#include "wasl/dg/scalar_advection/line_residual.hpp"
#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/history_csv.hpp"
#include "wasl/io/output_names.hpp"
#include "wasl/io/pvd_writer.hpp"
#include "wasl/io/pvtu_writer.hpp"
#include "wasl/io/vtu_writer.hpp"
#include "wasl/mesh/partition/partition.hpp"
#include "wasl/mesh/partition/rank_local_line_mesh.hpp"
#include "wasl/mesh/prepared/mesh_from_prepared.hpp"
#include "wasl/mesh/prepared/owned_mesh_slice.hpp"
#include "wasl/time/cfl_timestep.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {
namespace {

// Rank-local periodic 1D scalar advection on exactly two ranks.
//
// Each rank owns a contiguous block and imports its face-adjacent ghosts plus
// the periodic endpoint ghost across the domain: rank 0 imports the last
// element and rank 1 the first. Ghosts are exchanged at every SSP-RK3 stage.
// Ranks write owned VTU and checkpoint pieces; rank 0 writes the PVTU, PVD, and
// history.
class DistributedLineAdvection {
public:
  DistributedLineAdvection(const ScalarAdvectionCase &config, const CaseMesh &mesh,
                           const LineMesh &global)
      : config_(config), mesh_path_(mesh.path), prepared_(mesh.prepared), global_(global),
        rank_(world_rank()), rank_count_(world_size()),
        mode_count_(line_mode_count(config.polynomial_order)),
        last_element_(static_cast<std::int32_t>(global.num_elements() - 1)) {
    if (rank_count_ != 2)
      throw std::runtime_error(
          "simulation: MPI scalar-advection driver requires exactly two ranks");
    partition_ = partition_mesh(prepared_, rank_count_);
    assign_partition(prepared_, partition_);
    owned_ = rank_owned_elements(partition_, rank_);
    local_ = rank_local_line_mesh(global_, owned_, ghost_elements());
    coefficients_.assign(local_.mesh.num_elements() * mode_count_, 0.0);
  }

  SimulationSummary run() {
    if (config_.files.restart_checkpoint.empty())
      load_initial_state();
    else
      restore();
    const auto &run = config_.run;
    if (!std::isfinite(time_) || time_ < run.start_time || time_ > run.end_time ||
        step_ > run.max_steps)
      throw std::runtime_error(
          "simulation restart: checkpoint is outside the configured run range");

    std::filesystem::create_directories(config_.files.output_directory);
    owned_mesh_ = owned_mesh_slice(prepared_, rank_);
    history_.push_back({step_, time_, 0.0, global_mass()});
    write_snapshot();

    auto schedule = halo_schedule();
    HostStagedHaloExchange exchange(MPI_COMM_WORLD, schedule, mode_count_, 0, mode_count_,
                                    local_.global_to_local);
    const ScalarAdvection equation{config_.wavespeed};
    const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
      exchange.start(stage.data());
      exchange.finish(stage.data());
      const auto exterior = exterior_traces(stage);
      scalar_advection_residual(local_.mesh, equation, config_.polynomial_order, stage.data(),
                                exterior.data(), rhs.data());
    };
    const auto cfl_step =
        advection_cfl_timestep(global_, equation, config_.polynomial_order, config_.cfl);

    SspRk3Workspace workspace;
    while (time_ < run.end_time && step_ < run.max_steps) {
      const auto remaining = run.end_time - time_;
      const auto time_step = std::min(cfl_step > 0.0 ? cfl_step : remaining, remaining);
      if (time_step <= 0.0 || time_ + time_step == time_)
        throw std::runtime_error("simulation: timestep cannot advance time");
      ssp_rk3_step(time_step, residual, coefficients_, workspace);
      time_ += time_step;
      ++step_;
      const bool final = time_ >= run.end_time || step_ == run.max_steps;
      if (step_ % run.history_interval == 0 || final) {
        const auto mass = global_mass(); // collective on every rank
        if (rank_ == 0)
          history_.push_back({step_, time_, time_step, mass});
      }
      if (step_ % run.visualization_interval == 0 || final)
        write_snapshot();
      if (step_ % run.checkpoint_interval == 0 || final)
        write_restart();
    }
    if (time_ < run.end_time)
      throw std::runtime_error("simulation: max_steps reached before end_time");
    if (last_checkpoint_.empty())
      write_restart();
    world_barrier();
    if (rank_ == 0) {
      write_csv_history(config_.files.output_directory / "history.csv", history_);
      write_pvd(config_.files.output_directory / "solution.pvd", snapshots_);
    }
    return {step_, time_, last_checkpoint_};
  }

private:
  std::int32_t periodic_ghost() const { return rank_ == 0 ? last_element_ : 0; }
  std::int32_t periodic_send() const { return rank_ == 0 ? 0 : last_element_; }

  std::vector<std::int32_t> ghost_elements() const {
    auto ghosts = rank_ghost_elements(prepared_, partition_, rank_);
    if (std::find(ghosts.begin(), ghosts.end(), periodic_ghost()) == ghosts.end()) {
      ghosts.push_back(periodic_ghost());
      std::sort(ghosts.begin(), ghosts.end());
    }
    return ghosts;
  }

  HaloSchedule halo_schedule() const {
    auto schedule = build_halo_schedule(prepared_, partition_, rank_);
    const int peer = rank_ == 0 ? rank_count_ - 1 : 0;
    add_halo_send(schedule, peer, periodic_send());
    add_halo_recv(schedule, peer, periodic_ghost());
    return schedule;
  }

  std::size_t local_index(std::int32_t global) const {
    return static_cast<std::size_t>(local_.global_to_local.at(global));
  }

  // Owned elements occupy local slots 0..owned-1 in global order.
  void load_initial_state() {
    const auto initial = initial_coefficients(global_, config_);
    for (const auto global : owned_)
      std::copy_n(initial.data() + static_cast<std::size_t>(global) * mode_count_, mode_count_,
                  coefficients_.data() + local_index(global) * mode_count_);
  }

  void restore() {
    const auto checkpoint =
        read_checkpoint(rank_checkpoint_path(config_.files.restart_checkpoint, rank_));
    const std::vector<std::int64_t> expected(owned_.begin(), owned_.end());
    if (checkpoint.variable_count != 1 || checkpoint.mode_count != mode_count_ ||
        checkpoint.partition_rank != rank_ || checkpoint.partition_rank_count != rank_count_ ||
        checkpoint.global_element_count != global_.num_elements() ||
        checkpoint.element_indices != expected ||
        std::filesystem::weakly_canonical(checkpoint.prepared_mesh) != mesh_path_)
      throw std::runtime_error("simulation restart: checkpoint layout or mesh provenance mismatch");
    step_ = checkpoint.step;
    time_ = checkpoint.time;
    const auto restored = solver_coefficients(checkpoint);
    std::copy_n(restored.data(), owned_.size() * mode_count_, coefficients_.data());
  }

  double global_mass() const {
    double local_mass = 0.0;
    for (const auto global : owned_)
      local_mass += std::sqrt(2.0) *
                    element_map(global_, static_cast<std::size_t>(global)).jacobian() *
                    coefficients_[local_index(global) * mode_count_];
    return world_sum(local_mass);
  }

  // Boundary faces of the local mesh take the periodic neighbor's trace.
  std::vector<double> exterior_traces(const std::vector<double> &stage) const {
    const int order = config_.polynomial_order;
    const auto left_boundary = line_element_trace(stage, local_index(last_element_), order, 1.0);
    const auto right_boundary = line_element_trace(stage, local_index(0), order, -1.0);
    std::vector<double> exterior(local_.mesh.faces.size());
    for (std::size_t face = 0; face < local_.mesh.faces.size(); ++face) {
      const auto &value = local_.mesh.faces[face];
      if (value.is_boundary())
        exterior[face] = value.outward_normal < 0.0 ? left_boundary : right_boundary;
      else
        exterior[face] =
            line_element_trace(stage, static_cast<std::size_t>(value.right_element), order, -1.0);
    }
    return exterior;
  }

  void write_snapshot() {
    const auto base = numbered_step_name("solution", step_, "");
    std::vector<double> averages;
    averages.reserve(owned_.size());
    for (const auto global : owned_)
      averages.push_back(coefficients_[local_index(global) * mode_count_] / std::sqrt(2.0));
    write_vtu(config_.files.output_directory / rank_piece_name(base, rank_, ".vtu"), owned_mesh_,
              averages);
    world_barrier();
    if (rank_ != 0)
      return;
    std::vector<std::string> pieces;
    for (int source = 0; source < rank_count_; ++source)
      pieces.push_back(rank_piece_name(base, source, ".vtu"));
    const auto manifest = base + ".pvtu";
    write_pvtu(config_.files.output_directory / manifest, pieces);
    snapshots_.emplace_back(time_, manifest);
  }

  void write_restart() {
    const auto path = config_.files.output_directory /
                      rank_piece_name(numbered_step_name("checkpoint", step_, ""), rank_, ".h5");
    Checkpoint checkpoint;
    checkpoint.time = time_;
    checkpoint.step = step_;
    checkpoint.variable_count = 1;
    checkpoint.mode_count = mode_count_;
    checkpoint.element_count = owned_.size();
    checkpoint.prepared_mesh = mesh_path_.string();
    checkpoint.partition_rank = rank_;
    checkpoint.partition_rank_count = rank_count_;
    checkpoint.global_element_count = global_.num_elements();
    checkpoint.element_indices.assign(owned_.begin(), owned_.end());
    for (std::size_t mode = 0; mode < mode_count_; ++mode)
      for (const auto global : owned_)
        checkpoint.coefficients.push_back(coefficients_[local_index(global) * mode_count_ + mode]);
    write_checkpoint(path, checkpoint);
    checkpoints_.push_back(path);
    if (checkpoints_.size() > config_.run.keep_checkpoints) {
      std::filesystem::remove(checkpoints_.front());
      checkpoints_.pop_front();
    }
    last_checkpoint_ = path;
  }

  const ScalarAdvectionCase &config_;
  std::filesystem::path mesh_path_;
  PreparedMesh prepared_; // carries this run's two-rank partition
  const LineMesh &global_;
  int rank_;
  int rank_count_;
  std::size_t mode_count_;
  std::int32_t last_element_;
  Partition partition_;
  std::vector<std::int32_t> owned_;
  RankLocalLineMesh local_;
  PreparedMesh owned_mesh_;
  std::vector<double> coefficients_; // U[local element][mode]
  std::uint64_t step_ = 0;
  double time_ = config_.run.start_time;
  std::vector<HistoryRow> history_;
  std::vector<std::pair<double, std::string>> snapshots_;
  std::deque<std::filesystem::path> checkpoints_;
  std::filesystem::path last_checkpoint_;
};

} // namespace

SimulationSummary run_distributed_line_scalar_advection(const ScalarAdvectionCase &config,
                                                        const CaseMesh &mesh) {
  if (!all_periodic(config.boundaries))
    throw std::runtime_error("simulation: MPI scalar-advection boundaries must be periodic");
  const auto global = line_mesh_from_prepared(mesh.prepared);
  ensure_mpi_initialized();
  return DistributedLineAdvection(config, mesh, global).run();
}

} // namespace wasl

#else

namespace wasl {

SimulationSummary run_distributed_line_scalar_advection(const ScalarAdvectionCase &,
                                                        const CaseMesh &) {
  throw std::runtime_error("simulation: bridge was built without MPI support");
}

} // namespace wasl

#endif
