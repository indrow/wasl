#include "wasl/case/case_validation.hpp"

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace wasl {
namespace {

// Kernels are tabulated and tested through P3.
constexpr int maximum_polynomial_order = 3;

bool valid_files(const CaseFiles &files) {
  return !files.prepared_mesh.empty() && !files.output_directory.empty();
}

bool valid_polynomial_order(int order) { return order >= 0 && order <= maximum_polynomial_order; }

bool valid_run_control(const RunControl &run) {
  return std::isfinite(run.start_time) && std::isfinite(run.end_time) &&
         run.end_time > run.start_time && run.max_steps != 0 && run.history_interval != 0 &&
         run.visualization_interval != 0 && run.checkpoint_interval != 0 &&
         run.keep_checkpoints != 0;
}

bool positive_finite(double value) { return std::isfinite(value) && value > 0.0; }
bool nonnegative_finite(double value) { return std::isfinite(value) && value >= 0.0; }

bool valid_interval(double lower, double upper) {
  return std::isfinite(lower) && std::isfinite(upper) && lower < upper;
}

bool valid_circle(double center_x, double center_y, double radius) {
  return std::isfinite(center_x) && std::isfinite(center_y) && positive_finite(radius);
}

bool valid_primitive(const EulerPrimitive &value) {
  return positive_finite(value.density) && std::isfinite(value.velocity_x) &&
         std::isfinite(value.velocity_y) && positive_finite(value.pressure);
}

void validate_euler_region(const EulerRegion &region, int dimensions) {
  if (!valid_primitive(region.value))
    throw std::runtime_error("simulation: invalid Euler initial primitive state");
  if (region.shape == EulerRegion::Shape::interval &&
      (!valid_interval(region.lower, region.upper) || dimensions != 1))
    throw std::runtime_error("simulation: invalid Euler initial interval");
  if (region.shape == EulerRegion::Shape::circle &&
      (!valid_circle(region.center_x, region.center_y, region.radius) || dimensions != 2))
    throw std::runtime_error("simulation: invalid Euler initial circle");
}

// Names are unique and nonempty; an inflow value must be finite.
void validate_boundaries(const NamedBoundaries &boundaries) {
  std::set<std::string> names;
  for (const auto &boundary : boundaries) {
    if (boundary.name.empty() || !names.insert(boundary.name).second)
      throw std::runtime_error("simulation: boundary names must be unique and nonempty");
    if (!std::isfinite(boundary.condition.value))
      throw std::runtime_error("simulation: boundary " + boundary.name + " value must be finite");
  }
}

} // namespace

void validate_scalar_advection_case(const ScalarAdvectionCase &config) {
  if (!valid_files(config.files) || !valid_polynomial_order(config.polynomial_order) ||
      !std::isfinite(config.wavespeed) || !positive_finite(config.cfl) ||
      !std::isfinite(config.initial_default) || !std::isfinite(config.initial_sine_amplitude) ||
      !valid_run_control(config.run))
    throw std::runtime_error("simulation: invalid scalar-advection case configuration");
  if (config.dimensions < 1 || config.dimensions > 2 || !std::isfinite(config.velocity_y))
    throw std::runtime_error("simulation: invalid scalar-advection dimension or velocity");
  for (const auto &region : config.initial_intervals) {
    if (!valid_interval(region.lower, region.upper) || !std::isfinite(region.value))
      throw std::runtime_error("simulation: invalid initial interval");
  }
  for (const auto &region : config.initial_circles) {
    if (!valid_circle(region.center_x, region.center_y, region.radius) ||
        !std::isfinite(region.value))
      throw std::runtime_error("simulation: invalid initial circle");
  }
  validate_boundaries(config.boundaries);
}

void validate_euler_case(const EulerCase &config) {
  if (!valid_files(config.files) || !valid_polynomial_order(config.polynomial_order) ||
      !std::isfinite(config.heat_capacity_ratio) || config.heat_capacity_ratio <= 1.0 ||
      !nonnegative_finite(config.sensor_threshold) || !nonnegative_finite(config.density_floor) ||
      !nonnegative_finite(config.pressure_floor) || !positive_finite(config.cfl) ||
      !valid_run_control(config.run) || (config.dimensions != 1 && config.dimensions != 2) ||
      config.initial_regions.empty())
    throw std::runtime_error("simulation: invalid Euler case configuration");
  int default_count = 0;
  for (const auto &region : config.initial_regions) {
    validate_euler_region(region, config.dimensions);
    if (region.shape == EulerRegion::Shape::default_region)
      ++default_count;
  }
  if (default_count != 1)
    throw std::runtime_error("simulation: Euler initial state requires one default region");
  validate_boundaries(config.boundaries);
  for (const auto &boundary : config.boundaries) {
    if (boundary.condition.kind == BoundaryKind::inflow)
      throw std::runtime_error("simulation: boundary " + boundary.name +
                               ": Euler inflow is not implemented");
  }
}

} // namespace wasl
