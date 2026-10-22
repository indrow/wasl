// Linear scalar advection U_t + (a U)_x = S(U). For a canonical left-to-right
// face, the upwind flux is a U_left when a >= 0 and a U_right otherwise.

#pragma once

namespace wasl {

// Physical parameters for linear scalar advection. The wavespeed a is constant.
struct ScalarAdvection {
  double wavespeed;
  // Variable count is one for a scalar equation; exposed so equation-generic
  // code can size state without special cases.
  static constexpr int variable_count = 1;
};

// Constant physical velocity for two-dimensional linear scalar advection,
// U_t + div(v U) = 0.
struct ScalarAdvection2D {
  double velocity_x;
  double velocity_y;
  static constexpr int variable_count = 1;
};

// Physical flux F(U) = a U in the positive x direction.
inline double scalar_advection_flux(const ScalarAdvection &eq, double state) noexcept {
  return eq.wavespeed * state;
}

// Upwind numerical flux for linear advection across a face whose canonical
// normal points from the left state to the right state.
//
// The transported flux magnitude is a * U_upwind, where U_upwind is the state on
// the side the wind blows from: the left state when a > 0 and the right state
// when a < 0. The sign of the result carries the canonical-normal direction so
// callers can accumulate it consistently for both adjacent elements.
inline double upwind_flux(const ScalarAdvection &eq, double state_left,
                          double state_right) noexcept {
  if (eq.wavespeed >= 0.0) {
    return eq.wavespeed * state_left;
  }
  return eq.wavespeed * state_right;
}

// Upwind flux through a face with unit normal (normal_x, normal_y), directed
// from the left state to the right state.
inline double upwind_normal_flux(const ScalarAdvection2D &eq, double normal_x, double normal_y,
                                 double state_left, double state_right) noexcept {
  const double normal_speed = eq.velocity_x * normal_x + eq.velocity_y * normal_y;
  return normal_speed * (normal_speed >= 0.0 ? state_left : state_right);
}

} // namespace wasl
