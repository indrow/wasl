#include "wasl/case/serial_case_runner.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <string>
#include <utility>

#include "wasl/io/checkpoint_file.hpp"
#include "wasl/io/history_csv.hpp"
#include "wasl/io/output_names.hpp"
#include "wasl/io/pvd_writer.hpp"
#include "wasl/io/vtu_writer.hpp"

namespace wasl {
namespace {

using State = SerialCaseModel::State;

struct RunPosition {
  std::uint64_t step = 0;
  double time = 0.0;
};

// Load the restart checkpoint, rejecting one written for another layout or
// mesh.
State restore(const CaseFiles &files, const CaseMesh &mesh, const SerialCaseModel &model,
              RunPosition &position) {
  const auto checkpoint = read_checkpoint(files.restart_checkpoint);
  if (checkpoint.variable_count != model.variable_count ||
      checkpoint.mode_count != model.mode_count ||
      checkpoint.element_count != mesh.prepared.element_count() ||
      checkpoint.partition_rank_count != 1 ||
      std::filesystem::weakly_canonical(checkpoint.prepared_mesh) != mesh.path)
    throw std::runtime_error("simulation restart: checkpoint layout or mesh provenance mismatch");
  position.step = checkpoint.step;
  position.time = checkpoint.time;
  return model.from_checkpoint(checkpoint.coefficients);
}

void require_in_range(const RunControl &run, const RunPosition &position) {
  if (!std::isfinite(position.time) || position.time < run.start_time ||
      position.time > run.end_time || position.step > run.max_steps)
    throw std::runtime_error("simulation restart: checkpoint is outside the configured run range");
}

// Collects history and snapshots, and rotates checkpoints on disk.
class SerialOutput {
public:
  SerialOutput(const CaseFiles &files, const CaseMesh &mesh, const SerialCaseModel &model,
               std::uint64_t keep_checkpoints)
      : files_(files), mesh_(mesh), model_(model), keep_checkpoints_(keep_checkpoints) {
    std::filesystem::create_directories(files_.output_directory);
  }

  void record_history(const RunPosition &position, double time_step, const State &state) {
    history_.push_back({position.step, position.time, time_step, model_.mass(state)});
  }

  void write_snapshot(const RunPosition &position, const State &state) {
    const auto name = numbered_step_name("solution", position.step, ".vtu");
    write_vtu(files_.output_directory / name, mesh_.prepared, model_.cell_averages(state),
              model_.variable_count);
    snapshots_.emplace_back(position.time, name);
  }

  void write_restart(const RunPosition &position, const State &state) {
    const auto path =
        files_.output_directory / numbered_step_name("checkpoint", position.step, ".h5");
    Checkpoint checkpoint;
    checkpoint.time = position.time;
    checkpoint.step = position.step;
    checkpoint.variable_count = model_.variable_count;
    checkpoint.mode_count = model_.mode_count;
    checkpoint.element_count = mesh_.prepared.element_count();
    checkpoint.prepared_mesh = mesh_.path.string();
    checkpoint.coefficients = model_.to_checkpoint(state);
    write_checkpoint(path, checkpoint);
    checkpoints_.push_back(path);
    if (checkpoints_.size() > keep_checkpoints_) {
      std::filesystem::remove(checkpoints_.front());
      checkpoints_.pop_front();
    }
    last_checkpoint_ = path;
  }

  void finish() const {
    write_csv_history(files_.output_directory / "history.csv", history_);
    write_pvd(files_.output_directory / "solution.pvd", snapshots_);
  }

  const std::filesystem::path &last_checkpoint() const { return last_checkpoint_; }

private:
  const CaseFiles &files_;
  const CaseMesh &mesh_;
  const SerialCaseModel &model_;
  std::uint64_t keep_checkpoints_;
  std::vector<HistoryRow> history_;
  std::vector<std::pair<double, std::string>> snapshots_;
  std::deque<std::filesystem::path> checkpoints_;
  std::filesystem::path last_checkpoint_;
};

// The CFL step, clipped so the run lands exactly on end_time.
double next_time_step(const RunControl &run, double time, double stable_step) {
  const double remaining = run.end_time - time;
  const double time_step = std::min(stable_step > 0.0 ? stable_step : remaining, remaining);
  if (time_step <= 0.0 || time + time_step == time)
    throw std::runtime_error("simulation: timestep cannot advance time");
  return time_step;
}

} // namespace

SimulationSummary run_serial_case(const CaseFiles &files, const RunControl &run,
                                  const CaseMesh &mesh, const SerialCaseModel &model) {
  RunPosition position{0, run.start_time};
  State state = files.restart_checkpoint.empty() ? model.initial_state()
                                                 : restore(files, mesh, model, position);
  require_in_range(run, position);

  SerialOutput output(files, mesh, model, run.keep_checkpoints);
  output.record_history(position, 0.0, state);
  output.write_snapshot(position, state);

  SspRk3Workspace workspace;
  while (position.time < run.end_time && position.step < run.max_steps) {
    const double time_step = next_time_step(run, position.time, model.stable_time_step(state));
    model.advance(time_step, state, workspace);
    position.time += time_step;
    ++position.step;

    const bool final = position.time >= run.end_time || position.step == run.max_steps;
    if (position.step % run.history_interval == 0 || final)
      output.record_history(position, time_step, state);
    if (position.step % run.visualization_interval == 0 || final)
      output.write_snapshot(position, state);
    if (position.step % run.checkpoint_interval == 0 || final)
      output.write_restart(position, state);
  }
  if (position.time < run.end_time)
    throw std::runtime_error("simulation: max_steps reached before end_time");
  if (output.last_checkpoint().empty())
    output.write_restart(position, state);
  output.finish();
  return {position.step, position.time, output.last_checkpoint()};
}

} // namespace wasl
