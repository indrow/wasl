#include "case/euler_drivers.hpp"

#include "wasl/case/euler_diagnostics.hpp"
#include "wasl/limiting/euler_modal_limiter.hpp"
#include "wasl/mesh/line/line_mesh.hpp"
#include "wasl/mesh/triangle/triangle_mesh.hpp"

namespace wasl {

template <typename Mesh> void add_euler_state_operations(SerialCaseModel &model, const Mesh &mesh) {
  const auto element_count = mesh.num_elements();
  const auto mode_count = model.mode_count;
  model.from_checkpoint = [](const SerialCaseModel::State &stored) { return stored; };
  model.to_checkpoint = [](const SerialCaseModel::State &state) { return state; };
  model.mass = [&mesh](const SerialCaseModel::State &state) {
    return euler_density_mass(mesh, state);
  };
  model.cell_averages = [=](const SerialCaseModel::State &state) {
    return euler_averages(state, element_count, mode_count);
  };
}

template void add_euler_state_operations(SerialCaseModel &, const LineMesh &);
template void add_euler_state_operations(SerialCaseModel &, const TriangleMesh &);

StateFunction euler_stage_limiter(const EulerCase &config, const Euler &equation,
                                  const EulerLimiterStencil &stencil) {
  if (!config.limiting)
    return [](std::vector<double> &) {};
  return [&config, equation, &stencil](std::vector<double> &state) {
    limit_euler_modal(stencil, equation, config.sensor_threshold, config.density_floor,
                      config.pressure_floor, state);
  };
}

} // namespace wasl
