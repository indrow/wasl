// Local Lax-Friedrichs (Rusanov) flux for ideal-gas Euler.

#pragma once

#include <algorithm>

#include "wasl/equations/euler/flux.hpp"

namespace wasl {

inline EulerConserved rusanov_flux(const Euler &equation, const EulerConserved &left,
                                   const EulerConserved &right, double normal_x,
                                   double normal_y) noexcept {
  const EulerConserved left_flux = euler_normal_flux(equation, left, normal_x, normal_y);
  const EulerConserved right_flux = euler_normal_flux(equation, right, normal_x, normal_y);
  const double speed =
      std::max(euler_maximum_characteristic_speed(equation, left, normal_x, normal_y),
               euler_maximum_characteristic_speed(equation, right, normal_x, normal_y));
  return {0.5 * (left_flux.density + right_flux.density - speed * (right.density - left.density)),
          0.5 * (left_flux.momentum_x + right_flux.momentum_x -
                 speed * (right.momentum_x - left.momentum_x)),
          0.5 * (left_flux.momentum_y + right_flux.momentum_y -
                 speed * (right.momentum_y - left.momentum_y)),
          0.5 * (left_flux.total_energy + right_flux.total_energy -
                 speed * (right.total_energy - left.total_energy))};
}

} // namespace wasl
