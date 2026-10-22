// Third-order strong-stability-preserving Runge-Kutta (SSP-RK3).
//
// Residual operators return dc/dt with the mass inverse applied, so a
// forward-Euler stage is c + dt * rhs(c). SSP-RK3 composes three such stages in
// Shu-Osher form:
//   c1 = c                + dt * rhs(c)
//   c2 = (3/4) c          + (1/4) (c1 + dt * rhs(c1))
//   c  = (1/3) c          + (2/3) (c2 + dt * rhs(c2))
// The caller owns spatial discretization and boundary conditions through the
// residual callback; the stepper owns only time composition. This is the CPU
// reference path, so std::function and std::vector are acceptable here.

#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace wasl {

// Fill rhs from a stage. The stage is mutable so distributed callers can update
// ghost entries before evaluating the residual.
using ResidualFunction = std::function<void(std::vector<double> &, std::vector<double> &)>;

// Post-process each completed stage in place, for example by limiting.
using StateFunction = std::function<void(std::vector<double> &)>;

// Buffers reused across steps; steady-state stepping is allocation-free.
struct SspRk3Workspace {
  std::vector<double> stage;
  std::vector<double> rhs;
  std::vector<double> euler; // forward-Euler substate
  std::vector<double> exterior;
};

// One step with a stage post-processor applied after every stage.
inline void ssp_rk3_step(double dt, const ResidualFunction &residual,
                         const StateFunction &post_stage, std::vector<double> &coefficients,
                         SspRk3Workspace &workspace) {
  const std::size_t n = coefficients.size();
  workspace.stage.resize(n);
  workspace.rhs.resize(n);
  workspace.euler.resize(n);

  residual(coefficients, workspace.rhs);
  for (std::size_t i = 0; i < n; ++i)
    workspace.euler[i] = coefficients[i] + dt * workspace.rhs[i];
  post_stage(workspace.euler);

  residual(workspace.euler, workspace.rhs);
  for (std::size_t i = 0; i < n; ++i)
    workspace.stage[i] =
        0.75 * coefficients[i] + 0.25 * (workspace.euler[i] + dt * workspace.rhs[i]);
  post_stage(workspace.stage);

  residual(workspace.stage, workspace.rhs);
  for (std::size_t i = 0; i < n; ++i)
    coefficients[i] =
        (1.0 / 3.0) * coefficients[i] + (2.0 / 3.0) * (workspace.stage[i] + dt * workspace.rhs[i]);
  post_stage(coefficients);
}

// One step without stage post-processing.
inline void ssp_rk3_step(double dt, const ResidualFunction &residual,
                         std::vector<double> &coefficients, SspRk3Workspace &workspace) {
  ssp_rk3_step(dt, residual, [](std::vector<double> &) {}, coefficients, workspace);
}

} // namespace wasl
