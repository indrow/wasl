#include "wasl/case/euler_case.hpp"

#include "case/euler_drivers.hpp"
#include "wasl/case/case_validation.hpp"

namespace wasl {

SimulationSummary run_euler_case(const EulerCase &config) {
  validate_euler_case(config);
  const auto mesh = load_case_mesh(config.files.prepared_mesh, config.dimensions);
  return config.dimensions == 1 ? run_line_euler(config, mesh) : run_triangle_euler(config, mesh);
}

} // namespace wasl
