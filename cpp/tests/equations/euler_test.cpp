#include <cmath>

#include "wasl/boundaries/inflow.hpp"
#include "wasl/boundaries/outflow.hpp"
#include "wasl/boundaries/slip_wall.hpp"
#include "wasl/equations/euler/flux.hpp"
#include "wasl/flux/hllc.hpp"
#include "wasl/flux/rusanov.hpp"

namespace {

bool close(double actual, double expected, double tolerance = 1e-13) {
  return std::fabs(actual - expected) <= tolerance;
}

bool same(const wasl::EulerConserved &actual, const wasl::EulerConserved &expected,
          double tolerance = 1e-13) {
  return close(actual.density, expected.density, tolerance) &&
         close(actual.momentum_x, expected.momentum_x, tolerance) &&
         close(actual.momentum_y, expected.momentum_y, tolerance) &&
         close(actual.total_energy, expected.total_energy, tolerance);
}

bool check_conversion_and_flux() {
  const wasl::Euler equation{1.4};
  const wasl::EulerPrimitive initial{1.2, 2.0, -0.5, 1.0};
  const wasl::EulerConserved state = wasl::conserved(equation, initial);
  const wasl::EulerPrimitive recovered = wasl::primitive(equation, state);
  if (!close(recovered.density, initial.density) ||
      !close(recovered.velocity_x, initial.velocity_x) ||
      !close(recovered.velocity_y, initial.velocity_y) ||
      !close(recovered.pressure, initial.pressure) || !wasl::is_physical(equation, state))
    return false;

  const wasl::EulerConserved flux_x = wasl::euler_flux_x(equation, state);
  const wasl::EulerConserved flux_y = wasl::euler_flux_y(equation, state);
  return same(flux_x, {2.4, 5.8, -1.2, 12.1}) && same(flux_y, {-0.6, -1.2, 1.3, -3.025});
}

bool check_consistent_numerical_fluxes() {
  const wasl::Euler equation{1.4};
  const wasl::EulerConserved state =
      wasl::conserved(equation, wasl::EulerPrimitive{1.2, 2.0, -0.5, 1.0});
  const wasl::EulerConserved expected = wasl::euler_normal_flux(equation, state, 0.6, -0.8);
  return same(wasl::rusanov_flux(equation, state, state, 0.6, -0.8), expected) &&
         same(wasl::hllc_flux(equation, state, state, 0.6, -0.8), expected, 2e-13);
}

bool check_boundaries() {
  const wasl::Euler equation{1.4};
  const wasl::EulerConserved state =
      wasl::conserved(equation, wasl::EulerPrimitive{1.5, 2.0, -3.0, 0.8});
  const wasl::EulerConserved wall = wasl::slip_wall_state(state, 0.0, 1.0);
  return same(wasl::euler_inflow_state(state), state) &&
         same(wasl::euler_outflow_state(state), state) && close(wall.density, state.density) &&
         close(wall.momentum_x, state.momentum_x) && close(wall.momentum_y, -state.momentum_y) &&
         close(wall.total_energy, state.total_energy);
}

} // namespace

int main() {
  return check_conversion_and_flux() && check_consistent_numerical_fluxes() && check_boundaries()
             ? 0
             : 1;
}
