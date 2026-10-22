// Ideal-gas Euler case: configuration and entry point.

#pragma once

#include <vector>

#include "wasl/case/case_files.hpp"
#include "wasl/case/named_boundary.hpp"
#include "wasl/case/run_control.hpp"
#include "wasl/case/simulation_summary.hpp"
#include "wasl/equations/euler/euler.hpp"

namespace wasl {

// A primitive state over a region: the default region fills everything that
// no interval (1D) or circle (2D) contains.
struct EulerRegion {
  enum class Shape { default_region, interval, circle };
  Shape shape = Shape::default_region;
  double lower = 0.0;
  double upper = 0.0;
  double center_x = 0.0;
  double center_y = 0.0;
  double radius = 0.0;
  EulerPrimitive value{};
};

struct EulerCase {
  CaseFiles files;
  RunControl run;
  bool limiting = false;
  double sensor_threshold = 1e-3;
  double density_floor = 1e-10;
  double pressure_floor = 1e-10;
  int dimensions = 1;
  int polynomial_order = 0;
  double heat_capacity_ratio = 1.4;
  double cfl = 0.0;
  std::vector<EulerRegion> initial_regions;
  // Empty means periodic. Outflow is supported on 1D runs; inflow is not.
  NamedBoundaries boundaries;
};

// Euler with the HLLC flux on a prepared 1D or 2D mesh.
SimulationSummary run_euler_case(const EulerCase &config);

} // namespace wasl
