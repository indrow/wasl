// Linear scalar-advection case: configuration and entry point.

#pragma once

#include <vector>

#include "wasl/case/case_files.hpp"
#include "wasl/case/named_boundary.hpp"
#include "wasl/case/run_control.hpp"
#include "wasl/case/simulation_summary.hpp"

namespace wasl {

struct ScalarInterval {
  double lower;
  double upper;
  double value;
};

struct ScalarCircle {
  double center_x;
  double center_y;
  double radius;
  double value;
};

// The initial state is a sine wave when initial_sine is set; otherwise the
// first containing interval (1D) or circle (2D), else initial_default.
struct ScalarAdvectionCase {
  CaseFiles files;
  RunControl run;
  bool mpi = false;
  int dimensions = 1;
  int polynomial_order = 0;
  double wavespeed = 0.0; // x velocity
  double velocity_y = 0.0;
  double cfl = 0.0;
  double initial_default = 0.0;
  bool initial_sine = false;
  double initial_sine_amplitude = 0.0;
  std::vector<ScalarInterval> initial_intervals;
  std::vector<ScalarCircle> initial_circles;
  // Empty means periodic. Inflow and outflow are supported on serial 1D runs.
  NamedBoundaries boundaries;
};

// Scalar advection on a prepared 1D or 2D mesh. Periodic 1D cases may run on
// two MPI ranks when built with MPI support.
SimulationSummary run_scalar_advection_case(const ScalarAdvectionCase &config);

} // namespace wasl
