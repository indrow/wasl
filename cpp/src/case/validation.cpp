#include "wasl/case/validation.hpp"

#include <cmath>
#include <stdexcept>

namespace wasl {

void validate_scalar_advection_case(const ScalarAdvectionCase &config) {
  if (config.prepared_mesh.empty() || config.output_directory.empty() ||
      config.polynomial_order < 0 || config.polynomial_order > 3 ||
      !std::isfinite(config.wavespeed) || !std::isfinite(config.cfl) || config.cfl <= 0.0 ||
      !std::isfinite(config.start_time) || !std::isfinite(config.end_time) ||
      !std::isfinite(config.initial_default) || !std::isfinite(config.initial_sine_amplitude) ||
      config.end_time <= config.start_time || config.max_steps == 0 ||
      config.history_interval == 0 || config.visualization_interval == 0 ||
      config.checkpoint_interval == 0 || config.keep_checkpoints == 0) {
    throw std::runtime_error("simulation: invalid scalar-advection case configuration");
  }
  if (config.dimensions < 1 || config.dimensions > 2 || !std::isfinite(config.velocity_y)) {
    throw std::runtime_error("simulation: invalid scalar-advection dimension or velocity");
  }
  for (const auto &region : config.initial_intervals) {
    if (!std::isfinite(region.lower) || !std::isfinite(region.upper) ||
        !std::isfinite(region.value) || region.lower >= region.upper) {
      throw std::runtime_error("simulation: invalid initial interval");
    }
  }
  for (const auto &region : config.initial_circles) {
    if (!std::isfinite(region.center_x) || !std::isfinite(region.center_y) ||
        !std::isfinite(region.radius) || !std::isfinite(region.value) || region.radius <= 0.0) {
      throw std::runtime_error("simulation: invalid initial circle");
    }
  }
}

void validate_euler_case(const EulerCase &config) {
  if (config.prepared_mesh.empty() || config.output_directory.empty() ||
      config.polynomial_order < 0 || config.polynomial_order > 3 ||
      !std::isfinite(config.heat_capacity_ratio) || config.heat_capacity_ratio <= 1.0 ||
      !std::isfinite(config.sensor_threshold) || config.sensor_threshold < 0.0 ||
      !std::isfinite(config.density_floor) || config.density_floor < 0.0 ||
      !std::isfinite(config.pressure_floor) || config.pressure_floor < 0.0 ||
      !std::isfinite(config.cfl) || config.cfl <= 0.0 || !std::isfinite(config.start_time) ||
      !std::isfinite(config.end_time) || config.end_time <= config.start_time ||
      config.max_steps == 0 || config.history_interval == 0 || config.visualization_interval == 0 ||
      config.checkpoint_interval == 0 || config.keep_checkpoints == 0 ||
      (config.dimensions != 1 && config.dimensions != 2) || config.initial_regions.empty()) {
    throw std::runtime_error("simulation: invalid Euler case configuration");
  }
  int default_count = 0;
  for (const auto &region : config.initial_regions) {
    const auto &value = region.value;
    if (!std::isfinite(value.density) || !std::isfinite(value.velocity_x) ||
        !std::isfinite(value.velocity_y) || !std::isfinite(value.pressure) ||
        value.density <= 0.0 || value.pressure <= 0.0) {
      throw std::runtime_error("simulation: invalid Euler initial primitive state");
    }
    if (region.shape == EulerRegion::Shape::default_region) {
      ++default_count;
    } else if (region.shape == EulerRegion::Shape::interval &&
               (!std::isfinite(region.lower) || !std::isfinite(region.upper) ||
                region.lower >= region.upper || config.dimensions != 1)) {
      throw std::runtime_error("simulation: invalid Euler initial interval");
    } else if (region.shape == EulerRegion::Shape::circle &&
               (!std::isfinite(region.center_x) || !std::isfinite(region.center_y) ||
                !std::isfinite(region.radius) || region.radius <= 0.0 || config.dimensions != 2)) {
      throw std::runtime_error("simulation: invalid Euler initial circle");
    }
  }
  if (default_count != 1)
    throw std::runtime_error("simulation: Euler initial state requires one default region");
}

} // namespace wasl
