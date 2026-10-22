#include "wasl/case/scalar_advection_case.hpp"

#include "case/scalar_advection_drivers.hpp"
#include "wasl/case/case_validation.hpp"

namespace wasl {

SimulationSummary run_scalar_advection_case(const ScalarAdvectionCase &config) {
  validate_scalar_advection_case(config);
  const auto mesh = load_case_mesh(config.files.prepared_mesh, config.dimensions);
  if (config.dimensions == 2)
    return run_triangle_scalar_advection(config, mesh);
  if (config.mpi)
    return run_distributed_line_scalar_advection(config, mesh);
  return run_line_scalar_advection(config, mesh);
}

} // namespace wasl
