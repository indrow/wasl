// Mesh-specific Euler drivers behind run_euler_case.

#pragma once

#include "wasl/case/case_mesh.hpp"
#include "wasl/case/euler_case.hpp"
#include "wasl/case/serial_case_runner.hpp"
#include "wasl/case/simulation_summary.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/limiting/euler_limiter_stencil.hpp"

namespace wasl {

SimulationSummary run_line_euler(const EulerCase &config, const CaseMesh &mesh);
SimulationSummary run_triangle_euler(const EulerCase &config, const CaseMesh &mesh);

// Checkpoints already use the solver layout U[variable][mode][element].
template <typename Mesh> void add_euler_state_operations(SerialCaseModel &model, const Mesh &mesh);

// The configured limiter as a stage post-processor; a no-op when disabled.
StateFunction euler_stage_limiter(const EulerCase &config, const Euler &equation,
                                  const EulerLimiterStencil &stencil);

} // namespace wasl
