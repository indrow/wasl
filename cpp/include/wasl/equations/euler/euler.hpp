// Two-dimensional ideal-gas Euler state and equation parameters.

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

} // namespace wasl
