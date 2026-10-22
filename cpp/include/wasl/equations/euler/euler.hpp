// Two-dimensional ideal-gas Euler state and equation parameters.
//
// Conserved variables follow NUMERICAL_CONVENTIONS.md:
//   0 density, 1 x-momentum, 2 y-momentum, 3 total energy.

#pragma once

namespace wasl {

struct EulerConserved {
  double density;
  double momentum_x;
  double momentum_y;
  double total_energy;
};

struct EulerPrimitive {
  double density;
  double velocity_x;
  double velocity_y;
  double pressure;
};

struct Euler {
  double heat_capacity_ratio;
  static constexpr int variable_count = 4;
};

// Conserved variable by its layout index.
inline double euler_component(const EulerConserved &state, int variable) noexcept {
  switch (variable) {
  case 0:
    return state.density;
  case 1:
    return state.momentum_x;
  case 2:
    return state.momentum_y;
  default:
    return state.total_energy;
  }
}

inline double &euler_component(EulerConserved &state, int variable) noexcept {
  switch (variable) {
  case 0:
    return state.density;
  case 1:
    return state.momentum_x;
  case 2:
    return state.momentum_y;
  default:
    return state.total_energy;
  }
}

} // namespace wasl
