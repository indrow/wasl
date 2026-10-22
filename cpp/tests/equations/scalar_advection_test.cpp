// Tests for the linear scalar advection equation and upwind numerical flux.
//
// Checks: physical flux linearity, upwind selection for a > 0 and a < 0, the
// stationary (a = 0) case, and the P0 finite-volume reduction. The upwind flux
// for linear advection is exact, so these properties hold to machine precision.

#include <cmath>

#include "wasl/equations/scalar_advection/scalar_advection.hpp"

namespace {

constexpr double kTol = 1e-15;

bool approx(double a, double b, double tol = kTol) { return std::fabs(a - b) <= tol; }

// Physical flux F(U) = a U is linear in U.
bool check_physical_flux() {
  const wasl::ScalarAdvection eq{2.5};
  if (!approx(wasl::scalar_advection_flux(eq, 0.0), 0.0))
    return false;
  if (!approx(wasl::scalar_advection_flux(eq, 1.0), 2.5))
    return false;
  if (!approx(wasl::scalar_advection_flux(eq, -3.0), -7.5))
    return false;
  return true;
}

// Positive wavespeed selects the left (upwind) state; negative selects the right.
bool check_upwind_selection() {
  const wasl::ScalarAdvection forward{1.0};
  if (!approx(wasl::upwind_flux(forward, 7.0, 3.0), 7.0)) // a > 0: from the left
    return false;

  const wasl::ScalarAdvection backward{-1.0};
  if (!approx(wasl::upwind_flux(backward, 7.0, 3.0), -3.0)) // a < 0: from the right
    return false;
  return true;
}

// A zero wavespeed carries no information across the face.
bool check_stationary() {
  const wasl::ScalarAdvection still{0.0};
  if (!approx(wasl::upwind_flux(still, 7.0, 3.0), 0.0))
    return false;
  return true;
}

// Boundary a = 0 is classified as upwind-from-the-left (a >= 0), matching the
// conventional tie-break for a nonnegative wavespeed.
bool check_zero_is_upwind_left() {
  const wasl::ScalarAdvection still{0.0};
  // a == 0 falls in the >= 0 branch, so the left state is selected; both sides
  // give a*U = 0, so this only documents the tie-break, not a distinct value.
  return approx(wasl::upwind_flux(still, 7.0, 3.0), wasl::scalar_advection_flux(still, 7.0));
}

// P0 finite-volume reduction: with one mode per element, the upwind flux must
// reduce to the first-order upwind finite-volume flux a * U_upwind.
bool check_finite_volume_reduction() {
  const wasl::ScalarAdvection eq{1.0};
  // Left-upwind case.
  if (!approx(wasl::upwind_flux(eq, 2.0, 5.0), 2.0))
    return false;
  // Right-upwind case.
  const wasl::ScalarAdvection neg{-1.0};
  if (!approx(wasl::upwind_flux(neg, 2.0, 5.0), -5.0))
    return false;
  return true;
}

} // namespace

int main() {
  if (!check_physical_flux())
    return 1;
  if (!check_upwind_selection())
    return 1;
  if (!check_stationary())
    return 1;
  if (!check_zero_is_upwind_left())
    return 1;
  if (!check_finite_volume_reduction())
    return 1;
  return 0;
}
