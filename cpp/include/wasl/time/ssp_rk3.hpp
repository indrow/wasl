// SSP-RK3 time integration for the 1D scalar-advection DG residual.
//
// The residual operator returns dc_i/dt directly (the orthonormal mass inverse
// 1/J is already applied), so a forward-Euler update is c_i += dt * rhs_i.
// SSP-RK3 composes three such forward-Euler stages in the Shu-Osher form:
//   c1 = c                + dt * rhs(c)
//   c2 = (3/4) c          + (1/4) (c1 + dt * rhs(c1))
//   c  = (1/3) c          + (2/3) (c2 + dt * rhs(c2))
// Each stage needs the exterior trace state derived from that stage's
// coefficients, so the stepper takes a callback that recomputes the exterior
// state for whatever coefficients it is handed. This keeps boundary-condition
// ownership in the caller (periodic, inflow, outflow) while the stepper owns
// only time composition.
//
// State layout is U[element][mode], matching the residual operator. The stepper
// is header-only: it is the CPU reference path, so std::function and
// std::vector are acceptable here.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "wasl/basis/reference_element.hpp"
#include "wasl/dg/residual_operator.hpp"
#include "wasl/equations/scalar_advection/scalar_advection.hpp"
#include "wasl/mesh/line_geometry.hpp"
#include "wasl/mesh/line_mesh.hpp"

namespace wasl {

// Stable CFL timestep for 1D scalar advection on a line mesh. The DG advection
// stability scaling is (2p + 1): the spectral radius of the scaled operator.
// Returns dt = cfl * min_h / (|a| * (2p + 1)). Returns 0.0 when the wavespeed
// is zero (no stepping needed) or the mesh is empty.
inline double advection_cfl_timestep(const LineMesh &mesh, const ScalarAdvection &eq,
                                     int polynomial_order, double cfl) {
  if (mesh.elements.empty() || eq.wavespeed == 0.0) {
    return 0.0;
  }
  double min_length = mesh.elements.front().x_right - mesh.elements.front().x_left;
  for (const auto &element : mesh.elements) {
    min_length = std::min(min_length, element.x_right - element.x_left);
  }
  return cfl * min_length / (std::fabs(eq.wavespeed) * (2.0 * polynomial_order + 1.0));
}

// Recompute the exterior trace state for every face from the given coefficients.
// The caller supplies this; it encodes the boundary condition (periodic, inflow,
// outflow). The vector is indexed by face and read by scalar_advection_residual.
using ExteriorFunction = std::function<std::vector<double>(const std::vector<double> &)>;
// A residual callback may update ghost entries in the supplied stage before
// filling rhs. This lets distributed callers exchange every RK stage without
// duplicating the Shu-Osher composition.
using ResidualFunction = std::function<void(std::vector<double> &, std::vector<double> &)>;
using StateFunction = std::function<void(std::vector<double> &)>;

// Workspace buffers reused across the three RK stages to avoid per-step
// allocation. Size it once for the coefficient count and pass it to each step.
struct SspRk3Workspace {
  std::vector<double> stage;    // working coefficients (c1, then c2)
  std::vector<double> rhs;      // residual dc_i/dt
  std::vector<double> euler;    // stage + dt * rhs (the forward-Euler substate)
  std::vector<double> exterior; // exterior traces for the current stage
};

// Grow each buffer to at least n doubles. No reallocation happens after the
// first call for a given n, so steady-state stepping is allocation-free.
inline void reserve_workspace(SspRk3Workspace &ws, std::size_t n) {
  ws.stage.resize(n);
  ws.rhs.resize(n);
  ws.euler.resize(n);
}

// Apply one SSP-RK3 step using a caller-owned residual operation. The callback
// is invoked exactly once for each of the three stage states.
inline void ssp_rk3_step(double dt, const ResidualFunction &residual_fn,
                         std::vector<double> &coefficients, SspRk3Workspace &ws) {
  const std::size_t n = coefficients.size();
  reserve_workspace(ws, n);

  residual_fn(coefficients, ws.rhs);
  for (std::size_t i = 0; i < n; ++i) {
    ws.euler[i] = coefficients[i] + dt * ws.rhs[i];
  }

  residual_fn(ws.euler, ws.rhs);
  for (std::size_t i = 0; i < n; ++i) {
    ws.stage[i] = ws.euler[i] + dt * ws.rhs[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    ws.euler[i] = 0.75 * coefficients[i] + 0.25 * ws.stage[i];
  }

  residual_fn(ws.euler, ws.rhs);
  for (std::size_t i = 0; i < n; ++i) {
    ws.stage[i] = ws.euler[i] + dt * ws.rhs[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    coefficients[i] = (1.0 / 3.0) * coefficients[i] + (2.0 / 3.0) * ws.stage[i];
  }
}

inline void ssp_rk3_step(double dt, const ResidualFunction &residual_fn,
                         const StateFunction &state_fn, std::vector<double> &coefficients,
                         SspRk3Workspace &ws) {
  const std::size_t n = coefficients.size();
  reserve_workspace(ws, n);
  residual_fn(coefficients, ws.rhs);
  for (std::size_t i = 0; i < n; ++i)
    ws.euler[i] = coefficients[i] + dt * ws.rhs[i];
  state_fn(ws.euler);
  residual_fn(ws.euler, ws.rhs);
  for (std::size_t i = 0; i < n; ++i)
    ws.stage[i] = 0.75 * coefficients[i] + 0.25 * (ws.euler[i] + dt * ws.rhs[i]);
  state_fn(ws.stage);
  residual_fn(ws.stage, ws.rhs);
  for (std::size_t i = 0; i < n; ++i)
    coefficients[i] = (1.0 / 3.0) * coefficients[i] + (2.0 / 3.0) * (ws.stage[i] + dt * ws.rhs[i]);
  state_fn(coefficients);
}

// Apply one SSP-RK3 step in place to coefficients, using the fixed timestep dt.
// exterior_fn recomputes the face exterior state from a coefficient vector.
// ws is caller-owned scratch, reserved to coefficients.size().
inline void ssp_rk3_step(const LineMesh &mesh, const ScalarAdvection &eq, int polynomial_order,
                         double dt, const ExteriorFunction &exterior_fn,
                         std::vector<double> &coefficients, SspRk3Workspace &ws) {
  auto residual_of = [&](std::vector<double> &c, std::vector<double> &rhs) {
    ws.exterior = exterior_fn(c);
    scalar_advection_residual(mesh, eq, polynomial_order, c.data(), ws.exterior.data(), rhs.data());
  };
  ssp_rk3_step(dt, residual_of, coefficients, ws);
}

} // namespace wasl
