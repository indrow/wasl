// Cell-average-preserving positivity limiter for modal Euler states.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "wasl/equations/euler/euler.hpp"
#include "wasl/equations/euler/primitive.hpp"

namespace wasl {

inline EulerConserved evaluate_euler_modal(const double *coefficients, std::size_t mode_count,
                                           const double *basis) noexcept {
  EulerConserved state{0.0, 0.0, 0.0, 0.0};
  double *components[] = {&state.density, &state.momentum_x, &state.momentum_y,
                          &state.total_energy};
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      *components[variable] +=
          coefficients[static_cast<std::size_t>(variable) * mode_count + mode] * basis[mode];
  return state;
}

inline EulerConserved interpolate_euler_state(const EulerConserved &average,
                                              const EulerConserved &point, double theta) noexcept {
  return {average.density + theta * (point.density - average.density),
          average.momentum_x + theta * (point.momentum_x - average.momentum_x),
          average.momentum_y + theta * (point.momentum_y - average.momentum_y),
          average.total_energy + theta * (point.total_energy - average.total_energy)};
}

// basis_values is packed [point][mode]. Mode zero is phi_0 = 1/sqrt(2).
// Returns the total high-mode scale applied to coefficients.
inline double limit_euler_positivity(const Euler &equation, double *coefficients,
                                     std::size_t mode_count, const double *basis_values,
                                     std::size_t point_count, double density_floor,
                                     double pressure_floor) {
  if (mode_count == 0 || density_floor < 0.0 || pressure_floor < 0.0)
    throw std::invalid_argument("Euler positivity limiter: invalid modes or floors");

  const double inverse_sqrt2 = 1.0 / std::sqrt(2.0);
  const EulerConserved average{
      coefficients[0] * inverse_sqrt2, coefficients[mode_count] * inverse_sqrt2,
      coefficients[2 * mode_count] * inverse_sqrt2, coefficients[3 * mode_count] * inverse_sqrt2};
  if (average.density <= density_floor || primitive(equation, average).pressure <= pressure_floor)
    throw std::invalid_argument("Euler positivity limiter: cell average is below a floor");

  double density_scale = 1.0;
  for (std::size_t point = 0; point < point_count; ++point) {
    const EulerConserved state =
        evaluate_euler_modal(coefficients, mode_count, basis_values + point * mode_count);
    if (state.density < density_floor)
      density_scale = std::min(density_scale, (average.density - density_floor) /
                                                  (average.density - state.density));
  }
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    for (std::size_t mode = 1; mode < mode_count; ++mode)
      coefficients[static_cast<std::size_t>(variable) * mode_count + mode] *= density_scale;

  double pressure_scale = 1.0;
  for (std::size_t point = 0; point < point_count; ++point) {
    const EulerConserved state =
        evaluate_euler_modal(coefficients, mode_count, basis_values + point * mode_count);
    if (primitive(equation, state).pressure >= pressure_floor)
      continue;
    double lower = 0.0;
    double upper = 1.0;
    for (int iteration = 0; iteration < 64; ++iteration) {
      const double middle = 0.5 * (lower + upper);
      if (primitive(equation, interpolate_euler_state(average, state, middle)).pressure >=
          pressure_floor)
        lower = middle;
      else
        upper = middle;
    }
    pressure_scale = std::min(pressure_scale, lower);
  }
  for (int variable = 0; variable < Euler::variable_count; ++variable)
    for (std::size_t mode = 1; mode < mode_count; ++mode)
      coefficients[static_cast<std::size_t>(variable) * mode_count + mode] *= pressure_scale;
  return density_scale * pressure_scale;
}

} // namespace wasl
