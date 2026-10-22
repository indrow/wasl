#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "wasl/core/simulation.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

inline std::string numbered_step_name(const char *prefix, std::uint64_t step,
                                      const char *extension) {
  std::ostringstream name;
  name << prefix << '-' << std::setfill('0') << std::setw(6) << step << extension;
  return name.str();
}

// Shared serial case loop: restart restore, scheduled I/O, and SSP-RK3 advance hooks.
//
// `coefficients` are the initial solver state when restart_checkpoint is empty.
// On restart, the checkpoint coefficients are loaded, then `adapt_state` may convert
// layout (scalar solver layout) or leave them as-is (Euler checkpoint layout).
// `checkpoint_state` maps solver state into the persisted coefficient layout.
template <typename Mesh, typename Mass, typename Averages, typename AdaptState,
          typename CheckpointState, typename Initialize, typename Advance, typename Timestep,
          typename Synchronize>
SimulationSummary run_serial_case_loop(
    const std::filesystem::path &output_directory, const std::filesystem::path &prepared_path,
    const PreparedMesh &prepared, const Mesh &mesh, std::size_t mode_count,
    std::size_t variable_count, double start_time, double end_time, std::uint64_t max_steps,
    std::uint64_t history_interval, std::uint64_t visualization_interval,
    std::uint64_t checkpoint_interval, std::uint64_t keep_checkpoints,
    std::vector<double> coefficients, const std::filesystem::path &restart_checkpoint,
    const Mass &mass, const Averages &averages, const AdaptState &adapt_state,
    const CheckpointState &checkpoint_state, const Initialize &initialize, const Advance &advance,
    const Timestep &timestep, const Synchronize &synchronize) {
  const auto element_count = mesh.num_elements();
  std::uint64_t step = 0;
  double time = start_time;
  if (!restart_checkpoint.empty()) {
    const auto checkpoint = read_checkpoint(restart_checkpoint);
    if (checkpoint.variable_count != variable_count || checkpoint.mode_count != mode_count ||
        checkpoint.element_count != element_count || checkpoint.partition_rank_count != 1 ||
        std::filesystem::weakly_canonical(checkpoint.prepared_mesh) != prepared_path) {
      throw std::runtime_error("simulation restart: checkpoint layout or mesh provenance mismatch");
    }
    step = checkpoint.step;
    time = checkpoint.time;
    coefficients = adapt_state(checkpoint.coefficients);
  }
  if (!std::isfinite(time) || time < start_time || time > end_time || step > max_steps) {
    throw std::runtime_error("simulation restart: checkpoint is outside the configured run range");
  }
  initialize(coefficients);

  std::filesystem::create_directories(output_directory);
  std::vector<HistoryRow> history;
  std::vector<std::pair<double, std::string>> snapshots;
  std::deque<std::filesystem::path> checkpoints;
  SimulationSummary summary{step, time, {}};

  const auto write_visualization = [&] {
    const auto name = numbered_step_name("solution", step, ".vtu");
    write_vtu(output_directory / name, prepared, averages(coefficients), variable_count);
    snapshots.emplace_back(time, name);
  };
  const auto write_restart = [&] {
    const auto path = output_directory / numbered_step_name("checkpoint", step, ".h5");
    Checkpoint checkpoint;
    checkpoint.time = time;
    checkpoint.step = step;
    checkpoint.variable_count = variable_count;
    checkpoint.mode_count = mode_count;
    checkpoint.element_count = element_count;
    checkpoint.prepared_mesh = prepared_path.string();
    checkpoint.coefficients = checkpoint_state(coefficients);
    write_checkpoint(path, checkpoint);
    checkpoints.push_back(path);
    if (checkpoints.size() > keep_checkpoints) {
      std::filesystem::remove(checkpoints.front());
      checkpoints.pop_front();
    }
    summary.checkpoint = path;
  };

  history.push_back({step, time, 0.0, mass(coefficients)});
  write_visualization();
  SspRk3Workspace workspace;
  while (time < end_time && step < max_steps) {
    const auto cfl_step = timestep(coefficients);
    const auto time_step = std::min(cfl_step > 0.0 ? cfl_step : end_time - time, end_time - time);
    if (time_step <= 0.0 || time + time_step == time) {
      throw std::runtime_error("simulation: timestep cannot advance time");
    }
    advance(time_step, coefficients, workspace);
    time += time_step;
    ++step;
    const bool final = time >= end_time || step == max_steps;
    const bool write_history = step % history_interval == 0 || final;
    const bool write_solution = step % visualization_interval == 0 || final;
    const bool write_checkpoint = step % checkpoint_interval == 0 || final;
    if (write_history || write_solution || write_checkpoint) {
      synchronize(coefficients);
    }
    if (write_history) {
      history.push_back({step, time, time_step, mass(coefficients)});
    }
    if (write_solution) {
      write_visualization();
    }
    if (write_checkpoint) {
      write_restart();
    }
  }
  if (time < end_time) {
    throw std::runtime_error("simulation: max_steps reached before end_time");
  }
  if (summary.checkpoint.empty()) {
    write_restart();
  }
  write_csv_history(output_directory / "history.csv", history);
  write_pvd(output_directory / "solution.pvd", snapshots);
  summary.step = step;
  summary.time = time;
  return summary;
}

} // namespace wasl
