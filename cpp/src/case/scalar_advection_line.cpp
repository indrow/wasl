#include "case/scalar_advection_drivers.hpp"

#include "case/line_boundary_resolution.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/mesh/prepared/mesh_from_prepared.hpp"
#include "wasl/time/cfl_timestep.hpp"
#include "wasl/time/line_advection_step.hpp"

namespace wasl {

SimulationSummary run_line_scalar_advection(const ScalarAdvectionCase &config,
                                            const CaseMesh &mesh) {
  const auto line = line_mesh_from_prepared(mesh.prepared);
  const ScalarAdvection equation{config.wavespeed};
  const auto exterior =
      make_line_boundary_exterior(line, config.polynomial_order,
                                  resolve_line_boundaries(mesh.prepared, line, config.boundaries));
  const auto cfl_step = advection_cfl_timestep(line, equation, config.polynomial_order, config.cfl);

  SerialCaseModel model;
  model.variable_count = 1;
  model.mode_count = line_mode_count(config.polynomial_order);
  add_scalar_state_operations(model, line);
  model.initial_state = [&] { return initial_coefficients(line, config); };
  model.stable_time_step = [cfl_step](const SerialCaseModel::State &) { return cfl_step; };
  model.advance = [&](double time_step, SerialCaseModel::State &state, SspRk3Workspace &workspace) {
    ssp_rk3_step(line, equation, config.polynomial_order, time_step, exterior, state, workspace);
  };
  return run_serial_case(config.files, config.run, mesh, model);
}

} // namespace wasl
