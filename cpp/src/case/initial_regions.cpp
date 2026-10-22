#include "case/initial_regions.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace wasl {
namespace {

bool inside_circle(double x, double y, double center_x, double center_y, double radius) {
  const auto dx = x - center_x;
  const auto dy = y - center_y;
  return dx * dx + dy * dy <= radius * radius;
}

} // namespace

double scalar_initial_value(const ScalarAdvectionCase &config, double x) {
  if (config.initial_sine)
    return config.initial_sine_amplitude * std::sin(2.0 * std::numbers::pi * x);
  for (const auto &region : config.initial_intervals) {
    if (x >= region.lower && x <= region.upper)
      return region.value;
  }
  return config.initial_default;
}

double scalar_initial_value(const ScalarAdvectionCase &config, double x, double y) {
  for (const auto &region : config.initial_circles) {
    if (inside_circle(x, y, region.center_x, region.center_y, region.radius))
      return region.value;
  }
  return config.initial_default;
}

// The first containing interval or circle wins; the default fills the rest.
EulerPrimitive euler_initial_value(const EulerCase &config, double x, double y) {
  const EulerRegion *fallback = nullptr;
  for (const auto &region : config.initial_regions) {
    switch (region.shape) {
    case EulerRegion::Shape::default_region:
      fallback = &region;
      break;
    case EulerRegion::Shape::interval:
      if (x >= region.lower && x <= region.upper)
        return region.value;
      break;
    case EulerRegion::Shape::circle:
      if (inside_circle(x, y, region.center_x, region.center_y, region.radius))
        return region.value;
      break;
    }
  }
  if (fallback == nullptr)
    throw std::runtime_error("simulation: Euler initial state requires a default region");
  return fallback->value;
}

} // namespace wasl
