#include "case/euler_drivers.hpp"

#include "case/line_boundary_resolution.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/line_euler_exterior.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/dg/euler/line_residual.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/mesh/prepared/mesh_from_prepared.hpp"
#include "wasl/time/cfl_timestep.hpp"

namespace wasl {

SimulationSummary run_line_euler(const EulerCase &config, const CaseMesh &mesh) {
  const auto line = line_mesh_from_prepared(mesh.prepared);
  const Euler equation{config.heat_capacity_ratio};
  const auto mode_count = line_mode_count(config.polynomial_order);
  const auto conditions = resolve_line_boundaries(mesh.prepared, line, config.boundaries);
  const auto ends =
      conditions.left.kind == BoundaryKind::periodic ? LineEnds::periodic : LineEnds::open;
  const auto stencil = make_euler_limiter_stencil(line, config.polynomial_order, ends);
  const auto limit = euler_stage_limiter(config, equation, stencil);
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    const auto exterior = line_euler_exterior(line, config.polynomial_order, conditions, stage);
    euler_line_residual(line, equation, hllc_flux, config.polynomial_order, stage.data(),
                        exterior.data(), rhs.data());
  };

  SerialCaseModel model;
  model.variable_count = Euler::variable_count;
  model.mode_count = mode_count;
  add_euler_state_operations(model, line);
  model.initial_state = [&] {
    auto state = initial_coefficients(line, config, equation);
    limit(state);
    return state;
  };
  model.stable_time_step = [&](const SerialCaseModel::State &state) {
    return euler_cfl_timestep(line, equation, state, mode_count, config.polynomial_order,
                              config.cfl);
  };
  model.advance = [&](double time_step, SerialCaseModel::State &state, SspRk3Workspace &workspace) {
    ssp_rk3_step(time_step, residual, limit, state, workspace);
  };
  return run_serial_case(config.files, config.run, mesh, model);
}

} // namespace wasl
