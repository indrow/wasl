// Tests for the affine line physical/reference mapping.
//
// Checks: round-trip map consistency, constant Jacobian, endpoint mapping, and
// the reference interval bounds.

#include <cmath>

#include "wasl/mesh/line/line_geometry.hpp"

namespace {

constexpr double kTol = 1e-12;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

bool check_round_trip() {
  const wasl::LineMap map{0.17, 0.43};
  const double xs[] = {-1.0, -0.5, 0.0, 0.37, 1.0};
  for (double xi : xs) {
    const double x = map.reference_to_physical(xi);
    const double xi_back = map.physical_to_reference(x);
    if (!approx(xi_back, xi))
      return false;
  }
  return true;
}

bool check_endpoints_and_centre() {
  const wasl::LineMap map{0.17, 0.43};
  if (!approx(map.reference_to_physical(-1.0), 0.17))
    return false;
  if (!approx(map.reference_to_physical(1.0), 0.43))
    return false;
  if (!approx(map.centre(), 0.30))
    return false;
  return true;
}

bool check_constant_jacobian() {
  // Irregular element: Jacobian is half the length and independent of xi.
  const wasl::LineMap map{0.71, 1.0};
  if (!approx(map.jacobian(), 0.145))
    return false;
  if (!approx(map.length(), 0.29))
    return false;
  return true;
}

// A constant field mapped from reference to physical stays constant: the map is
// affine, so integrating a constant over physical space equals length times the
// reference average.
bool check_constant_metric() {
  const wasl::LineMap map{0.0, 1.0};
  // integral of 1 over [0,1] in physical space = 1; reference integral is 2;
  // metric factor length/2 reconciles them: 2 * (1/2) = 1.
  const double reference_integral = 2.0; // integral of 1 over [-1,1]
  if (!approx(reference_integral * map.jacobian(), map.length()))
    return false;
  return true;
}

} // namespace

int main() {
  if (!check_round_trip())
    return 1;
  if (!check_endpoints_and_centre())
    return 1;
  if (!check_constant_jacobian())
    return 1;
  if (!check_constant_metric())
    return 1;
  return 0;
}
