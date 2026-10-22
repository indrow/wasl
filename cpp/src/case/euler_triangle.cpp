#include "case/euler_drivers.hpp"

#include <stdexcept>

#include "case/line_boundary_resolution.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/periodic_triangle_trace.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/dg/euler/triangle_residual.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/mesh/prepared/mesh_from_prepared.hpp"
#include "wasl/time/cfl_timestep.hpp"

namespace wasl {

SimulationSummary run_triangle_euler(const EulerCase &config, const CaseMesh &mesh) {
  if (!all_periodic(config.boundaries))
    throw std::runtime_error("simulation: 2D Euler boundaries must be periodic");
  const auto triangles = triangle_mesh_from_prepared(mesh.prepared);
  const Euler equation{config.heat_capacity_ratio};
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto trace = periodic_triangle_trace(triangles, config.polynomial_order);
  const auto stencil = make_euler_limiter_stencil(triangles, config.polynomial_order,
                                                  trace.elements, trace.point_count);
  const auto limit = euler_stage_limiter(config, equation, stencil);
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    const auto exterior = periodic_triangle_euler_exterior(trace, stage, triangles.num_elements(),
                                                           mode_count, triangles.faces.size());
    euler_triangle_residual(triangles, equation, hllc_flux, config.polynomial_order, stage.data(),
                            exterior.data(), rhs.data());
  };

  SerialCaseModel model;
  model.variable_count = Euler::variable_count;
  model.mode_count = mode_count;
  add_euler_state_operations(model, triangles);
  model.initial_state = [&] {
    auto state = initial_coefficients(triangles, config, equation);
    limit(state);
    return state;
  };
  model.stable_time_step = [&](const SerialCaseModel::State &state) {
    return euler_cfl_timestep(triangles, equation, state, mode_count, config.polynomial_order,
                              config.cfl);
  };
  model.advance = [&](double time_step, SerialCaseModel::State &state, SspRk3Workspace &workspace) {
    ssp_rk3_step(time_step, residual, limit, state, workspace);
  };
  return run_serial_case(config.files, config.run, mesh, model);
}

} // namespace wasl
