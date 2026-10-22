// HLLC approximate Riemann flux for ideal-gas Euler in an arbitrary direction.

#pragma once

#include <algorithm>
#include <cmath>

#include "wasl/equations/euler/flux.hpp"

namespace wasl {
namespace detail {

inline EulerConserved hllc_star_state(const EulerConserved &state,
                                      const EulerPrimitive &primitive_state, double normal_x,
                                      double normal_y, double wave_speed,
                                      double contact_speed) noexcept {
  const double normal_velocity =
      primitive_state.velocity_x * normal_x + primitive_state.velocity_y * normal_y;
  const double star_density =
      state.density * (wave_speed - normal_velocity) / (wave_speed - contact_speed);
  const double velocity_shift = contact_speed - normal_velocity;
  const double specific_star_energy =
      state.total_energy / state.density +
      velocity_shift * (contact_speed + primitive_state.pressure /
                                            (state.density * (wave_speed - normal_velocity)));
  return {star_density, star_density * (primitive_state.velocity_x + velocity_shift * normal_x),
          star_density * (primitive_state.velocity_y + velocity_shift * normal_y),
          star_density * specific_star_energy};
}

inline EulerConserved add_scaled_jump(const EulerConserved &flux, double speed,
                                      const EulerConserved &star,
                                      const EulerConserved &state) noexcept {
  return {flux.density + speed * (star.density - state.density),
          flux.momentum_x + speed * (star.momentum_x - state.momentum_x),
          flux.momentum_y + speed * (star.momentum_y - state.momentum_y),
          flux.total_energy + speed * (star.total_energy - state.total_energy)};
}

} // namespace detail

inline EulerConserved hllc_flux(const Euler &equation, const EulerConserved &left,
                                const EulerConserved &right, double normal_x,
                                double normal_y) noexcept {
  const EulerPrimitive left_value = primitive(equation, left);
  const EulerPrimitive right_value = primitive(equation, right);
  const double left_normal_velocity =
      left_value.velocity_x * normal_x + left_value.velocity_y * normal_y;
  const double right_normal_velocity =
      right_value.velocity_x * normal_x + right_value.velocity_y * normal_y;
  const double left_sound_speed =
      std::sqrt(equation.heat_capacity_ratio * left_value.pressure / left_value.density);
  const double right_sound_speed =
      std::sqrt(equation.heat_capacity_ratio * right_value.pressure / right_value.density);
  const double left_wave =
      std::min(left_normal_velocity - left_sound_speed, right_normal_velocity - right_sound_speed);
  const double right_wave =
      std::max(left_normal_velocity + left_sound_speed, right_normal_velocity + right_sound_speed);
  const EulerConserved left_flux = euler_normal_flux(equation, left, normal_x, normal_y);
  if (left_wave >= 0.0)
    return left_flux;
  const EulerConserved right_flux = euler_normal_flux(equation, right, normal_x, normal_y);
  if (right_wave <= 0.0)
    return right_flux;

  const double contact_speed =
      (right_value.pressure - left_value.pressure +
       left.density * left_normal_velocity * (left_wave - left_normal_velocity) -
       right.density * right_normal_velocity * (right_wave - right_normal_velocity)) /
      (left.density * (left_wave - left_normal_velocity) -
       right.density * (right_wave - right_normal_velocity));
  if (contact_speed >= 0.0) {
    const EulerConserved star =
        detail::hllc_star_state(left, left_value, normal_x, normal_y, left_wave, contact_speed);
    return detail::add_scaled_jump(left_flux, left_wave, star, left);
  }
  const EulerConserved star =
      detail::hllc_star_state(right, right_value, normal_x, normal_y, right_wave, contact_speed);
  return detail::add_scaled_jump(right_flux, right_wave, star, right);
}

} // namespace wasl
