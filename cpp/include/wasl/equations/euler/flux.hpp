// Ideal-gas Euler physical fluxes and characteristic speed.

#pragma once

#include <cmath>

#include "wasl/equations/euler/primitive.hpp"

namespace wasl {

inline EulerConserved euler_flux_x(const Euler &equation, const EulerConserved &state) noexcept {
  const EulerPrimitive value = primitive(equation, state);
  return {state.momentum_x, state.momentum_x * value.velocity_x + value.pressure,
          state.momentum_y * value.velocity_x,
          (state.total_energy + value.pressure) * value.velocity_x};
}

inline EulerConserved euler_flux_y(const Euler &equation, const EulerConserved &state) noexcept {
  const EulerPrimitive value = primitive(equation, state);
  return {state.momentum_y, state.momentum_x * value.velocity_y,
          state.momentum_y * value.velocity_y + value.pressure,
          (state.total_energy + value.pressure) * value.velocity_y};
}

inline EulerConserved euler_normal_flux(const Euler &equation, const EulerConserved &state,
                                        double normal_x, double normal_y) noexcept {
  const EulerConserved flux_x = euler_flux_x(equation, state);
  const EulerConserved flux_y = euler_flux_y(equation, state);
  return {normal_x * flux_x.density + normal_y * flux_y.density,
          normal_x * flux_x.momentum_x + normal_y * flux_y.momentum_x,
          normal_x * flux_x.momentum_y + normal_y * flux_y.momentum_y,
          normal_x * flux_x.total_energy + normal_y * flux_y.total_energy};
}

inline double euler_maximum_characteristic_speed(const Euler &equation, const EulerConserved &state,
                                                 double normal_x, double normal_y) noexcept {
  const EulerPrimitive value = primitive(equation, state);
  return std::fabs(normal_velocity(value, normal_x, normal_y)) + sound_speed(equation, value);
}

} // namespace wasl
