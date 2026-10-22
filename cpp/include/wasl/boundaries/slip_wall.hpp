// Inviscid slip-wall exterior state for Euler.

#pragma once

#include "wasl/equations/euler/euler.hpp"

namespace wasl {

inline EulerConserved slip_wall_state(const EulerConserved &interior, double normal_x,
                                      double normal_y) noexcept {
  const double normal_momentum = interior.momentum_x * normal_x + interior.momentum_y * normal_y;
  return {interior.density, interior.momentum_x - 2.0 * normal_momentum * normal_x,
          interior.momentum_y - 2.0 * normal_momentum * normal_y, interior.total_energy};
}

} // namespace wasl
