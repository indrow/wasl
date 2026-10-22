// Mesh- and transport-specific scalar-advection drivers behind
// run_scalar_advection_case.

#pragma once

#include "wasl/case/case_mesh.hpp"
#include "wasl/case/scalar_advection_case.hpp"
#include "wasl/case/serial_case_runner.hpp"
#include "wasl/case/simulation_summary.hpp"

namespace wasl {

SimulationSummary run_line_scalar_advection(const ScalarAdvectionCase &config,
                                            const CaseMesh &mesh);
SimulationSummary run_triangle_scalar_advection(const ScalarAdvectionCase &config,
                                                const CaseMesh &mesh);

// Two-rank MPI driver; throws when built without MPI support.
SimulationSummary run_distributed_line_scalar_advection(const ScalarAdvectionCase &config,
                                                        const CaseMesh &mesh);

// Layout conversions and diagnostics shared by both serial scalar drivers.
template <typename Mesh> void add_scalar_state_operations(SerialCaseModel &model, const Mesh &mesh);

} // namespace wasl
