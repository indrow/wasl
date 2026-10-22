// The single-rank case loop: restart, CFL stepping, and scheduled output.
//
// Each case type describes itself as a SerialCaseModel. The runner owns time,
// step counting, history, snapshots, and checkpoint rotation; the model owns
// the state layout and the numerical advance. Runtime polymorphism stays at
// this orchestration level (ARCHITECTURE.md).

#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "wasl/case/case_files.hpp"
#include "wasl/case/case_mesh.hpp"
#include "wasl/case/run_control.hpp"
#include "wasl/case/simulation_summary.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

struct SerialCaseModel {
  using State = std::vector<double>;

  std::size_t variable_count = 0;
  std::size_t mode_count = 0;
  std::function<State()> initial_state;                  // fresh starts only
  std::function<State(const State &)> from_checkpoint;   // persisted -> solver layout
  std::function<State(const State &)> to_checkpoint;     // solver -> persisted layout
  std::function<double(const State &)> mass;             // history diagnostic
  std::function<State(const State &)> cell_averages;     // packed [element][variable]
  std::function<double(const State &)> stable_time_step; // <= 0 means unrestricted
  std::function<void(double, State &, SspRk3Workspace &)> advance;
};

SimulationSummary run_serial_case(const CaseFiles &files, const RunControl &run,
                                  const CaseMesh &mesh, const SerialCaseModel &model);

} // namespace wasl
