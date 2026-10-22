#include "case/scalar_advection_drivers.hpp"

#include <stdexcept>

#include "case/line_boundary_resolution.hpp"
#include "wasl/basis/mode_count.hpp"
#include "wasl/boundaries/periodic_triangle_trace.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/dg/scalar_advection/triangle_residual.hpp"
#include "wasl/mesh/prepared/mesh_from_prepared.hpp"
#include "wasl/time/cfl_timestep.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {

SimulationSummary run_triangle_scalar_advection(const ScalarAdvectionCase &config,
                                                const CaseMesh &mesh) {
  if (config.mpi)
    throw std::runtime_error("simulation: 2D case-level MPI is not implemented");
  if (!all_periodic(config.boundaries))
    throw std::runtime_error("simulation: 2D scalar-advection boundaries must be periodic");
  const auto triangles = triangle_mesh_from_prepared(mesh.prepared);
  const ScalarAdvection2D equation{config.wavespeed, config.velocity_y};
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto trace = periodic_triangle_trace(triangles, config.polynomial_order);
  const auto cfl_step = triangle_scalar_cfl_timestep(triangles, config.wavespeed, config.velocity_y,
                                                     config.polynomial_order, config.cfl);
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    const auto exterior = periodic_triangle_exterior(trace, stage, mode_count);
    triangle_scalar_advection_residual(triangles, equation, config.polynomial_order, stage.data(),
                                       exterior.data(), rhs.data());
  };

  SerialCaseModel model;
  model.variable_count = 1;
  model.mode_count = mode_count;
  add_scalar_state_operations(model, triangles);
  model.initial_state = [&] { return initial_coefficients(triangles, config); };
  model.stable_time_step = [cfl_step](const SerialCaseModel::State &) { return cfl_step; };
  model.advance = [&](double time_step, SerialCaseModel::State &state, SspRk3Workspace &workspace) {
    ssp_rk3_step(time_step, residual, state, workspace);
  };
  return run_serial_case(config.files, config.run, mesh, model);
}

} // namespace wasl
