// Allocation-free ideal-gas conserved and primitive state conversions.

#pragma once

#include <cmath>

#include "wasl/equations/euler/euler.hpp"

namespace wasl {

inline EulerConserved conserved(const Euler &equation, const EulerPrimitive &state) noexcept {
  const double kinetic_energy =
      0.5 * state.density *
      (state.velocity_x * state.velocity_x + state.velocity_y * state.velocity_y);
  return {state.density, state.density * state.velocity_x, state.density * state.velocity_y,
          state.pressure / (equation.heat_capacity_ratio - 1.0) + kinetic_energy};
}

inline EulerPrimitive primitive(const Euler &equation, const EulerConserved &state) noexcept {
  const double inverse_density = 1.0 / state.density;
  const double velocity_x = state.momentum_x * inverse_density;
  const double velocity_y = state.momentum_y * inverse_density;
  const double kinetic_energy =
      0.5 * (state.momentum_x * velocity_x + state.momentum_y * velocity_y);
  return {state.density, velocity_x, velocity_y,
          (equation.heat_capacity_ratio - 1.0) * (state.total_energy - kinetic_energy)};
}

inline double sound_speed(const Euler &equation, const EulerPrimitive &state) noexcept {
  return std::sqrt(equation.heat_capacity_ratio * state.pressure / state.density);
}

inline double normal_velocity(const EulerPrimitive &state, double normal_x,
                              double normal_y) noexcept {
  return state.velocity_x * normal_x + state.velocity_y * normal_y;
}

inline bool is_physical(const Euler &equation, const EulerConserved &state) noexcept {
  return state.density > 0.0 && primitive(equation, state).pressure > 0.0;
}

} // namespace wasl
