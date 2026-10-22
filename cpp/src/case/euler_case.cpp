#include "wasl/core/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "wasl/basis/reference_element.hpp"
#include "wasl/case/coefficient_layout.hpp"
#include "wasl/case/initial_state.hpp"
#include "wasl/case/periodic_boundary.hpp"
#include "wasl/case/prepared_mesh.hpp"
#include "wasl/case/run_schedule.hpp"
#include "wasl/case/validation.hpp"
#include "wasl/dg/euler_line_residual.hpp"
#include "wasl/dg/euler_triangle_residual.hpp"
#include "wasl/equations/euler/euler.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/io/case_io.hpp"
#include "wasl/limiting/euler_modal_limiter.hpp"
#include "wasl/mesh/line_mesh.hpp"
#include "wasl/mesh/triangle_mesh.hpp"
#include "wasl/time/ssp_rk3.hpp"

namespace wasl {
namespace {

template <typename Mesh, typename Initialize, typename Advance, typename Timestep,
          typename Synchronize>
SimulationSummary run_serial_euler_case(const EulerCase &config,
                                        const std::filesystem::path &prepared_path,
                                        const PreparedMesh &prepared, const Mesh &mesh,
                                        std::size_t mode_count, std::vector<double> initial,
                                        const Initialize &initialize, const Advance &advance,
                                        const Timestep &timestep, const Synchronize &synchronize) {
  const auto element_count = mesh.num_elements();
  const auto mass = [&](const std::vector<double> &state) {
    return euler_density_mass(mesh, state, mode_count);
  };
  const auto averages = [&](const std::vector<double> &state) {
    return euler_averages(state, element_count, mode_count);
  };
  const auto adapt_state = [](const std::vector<double> &stored) { return stored; };
  const auto checkpoint_state = [](const std::vector<double> &state) { return state; };
  return run_serial_case_loop(config.output_directory, prepared_path, prepared, mesh, mode_count,
                              static_cast<std::size_t>(Euler::variable_count), config.start_time,
                              config.end_time, config.max_steps, config.history_interval,
                              config.visualization_interval, config.checkpoint_interval,
                              config.keep_checkpoints, std::move(initial),
                              config.restart_checkpoint, mass, averages, adapt_state,
                              checkpoint_state, initialize, advance, timestep, synchronize);
}

} // namespace

SimulationSummary run_euler_case(const EulerCase &config) {
  validate_euler_case(config);
  const auto prepared_path = std::filesystem::weakly_canonical(config.prepared_mesh);
  const auto prepared = read_prepared_mesh(prepared_path);
  if (config.dimensions != prepared.dimension)
    throw std::runtime_error("simulation: case dimension does not match prepared mesh");
  const Euler equation{config.heat_capacity_ratio};
  if (config.dimensions == 1) {
    const auto mesh = line_mesh_from_prepared(prepared);
    const auto mode_count = line_mode_count(config.polynomial_order);
    const auto limiter = make_euler_limiter_data(mesh, config.polynomial_order);
    auto initial = initial_coefficients(mesh, config, equation);
    if (config.limiting)
      limit_euler_modal(limiter, equation, config.sensor_threshold, config.density_floor,
                        config.pressure_floor, initial);
    const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
      const auto exterior = periodic_line_euler_exterior(mesh, config.polynomial_order, stage);
      euler_line_residual(mesh, equation, hllc_flux, config.polynomial_order, stage.data(),
                          exterior.data(), rhs.data());
    };
    const auto timestep = [&](const std::vector<double> &stage) {
      return euler_cfl_timestep(mesh, equation, stage, mode_count, config.polynomial_order,
                                config.cfl);
    };
    const auto initialize = [](const std::vector<double> &) {};
    const StateFunction limit = [&](std::vector<double> &stage) {
      if (config.limiting)
        limit_euler_modal(limiter, equation, config.sensor_threshold, config.density_floor,
                          config.pressure_floor, stage);
    };
    const auto advance = [&](double time_step, std::vector<double> &stage,
                             SspRk3Workspace &workspace) {
      ssp_rk3_step(time_step, residual, limit, stage, workspace);
    };
    const auto synchronize = [](std::vector<double> &) {};
    return run_serial_euler_case(config, prepared_path, prepared, mesh, mode_count,
                                 std::move(initial), initialize, advance, timestep, synchronize);
  }

  const auto mesh = triangle_mesh_from_prepared(prepared);
  const auto mode_count = triangle_mode_count(config.polynomial_order);
  const auto trace = periodic_triangle_trace(mesh, config.polynomial_order);
  const auto limiter =
      make_euler_limiter_data(mesh, config.polynomial_order, trace.elements, trace.point_count);
  auto initial = initial_coefficients(mesh, config, equation);
  if (config.limiting)
    limit_euler_modal(limiter, equation, config.sensor_threshold, config.density_floor,
                      config.pressure_floor, initial);
  const ResidualFunction residual = [&](std::vector<double> &stage, std::vector<double> &rhs) {
    const auto exterior = periodic_triangle_euler_exterior(trace, stage, mesh.num_elements(),
                                                           mode_count, mesh.faces.size());
    euler_triangle_residual(mesh, equation, hllc_flux, config.polynomial_order, stage.data(),
                            exterior.data(), rhs.data());
  };
  const auto timestep = [&](const std::vector<double> &stage) {
    return euler_cfl_timestep(mesh, equation, stage, mode_count, config.polynomial_order,
                              config.cfl);
  };
  const auto initialize = [](const std::vector<double> &) {};
  const StateFunction limit = [&](std::vector<double> &stage) {
    if (config.limiting)
      limit_euler_modal(limiter, equation, config.sensor_threshold, config.density_floor,
                        config.pressure_floor, stage);
  };
  const auto advance = [&](double time_step, std::vector<double> &stage,
                           SspRk3Workspace &workspace) {
    ssp_rk3_step(time_step, residual, limit, stage, workspace);
  };
  const auto synchronize = [](std::vector<double> &) {};
  return run_serial_euler_case(config, prepared_path, prepared, mesh, mode_count,
                               std::move(initial), initialize, advance, timestep, synchronize);
}

} // namespace wasl
